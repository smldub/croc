use r_libcint as cint;
use croc_parse::{parse_nwchem, BASIS_MAP};
use statrs::function::gamma::gamma;
use itertools::{Itertools,iproduct};
use crate::ATOMIC_NUMBER_MAP;
use ndarray::{array,Array,Array2,Array4, s, Order};

const BASIS_EXCEPTIONS: [char;2] = ['*','+'];

// Might need a built Mol object that no longer allow modification to make sure the results match
// the stuff
pub struct Mol {
    /// The jth atom's name and coords is stored as 
    /// atoms[j] = (name,(x,y,z))
    atoms: Vec<(String,f64,f64,f64)>,
    /// Name of basis (see TODO for set of available bases)
    basis: String,
    /// c_dat contains the correctly formatted data for libcint
    c_dat: cint::CINTR2CDATA,
}

impl Mol {
    pub fn new() -> Self {
        Self{
            atoms: Vec::new(),
            basis: String::new(),
            c_dat: cint::CINTR2CDATA::new(),
        }
    }
    pub fn set_basis(&mut self, basis_name: &str) -> Result<(),&str> {
        let temp = basis_name.trim()
                             .to_lowercase()
                             .chars()
                             .filter(|x| x.is_alphanumeric() | &BASIS_EXCEPTIONS.contains(x))
                             .collect::<String>();
        if BASIS_MAP.contains_key(temp.as_str()) {
            self.basis = temp;
            Ok(())
        }
        else {
            Err("Basis is not included in std library")
        }
    }
    pub fn set_atoms_str(&mut self, atm_coords: &str) -> Result<(),&str>{
        /// Format = "a1 x1 y1 z1; a2 x2 y2 z2;...."
        let input = atm_coords.trim().split_whitespace().collect::<Vec<&str>>();
        self.atoms = (0..input.len()/4).filter(|i| ATOMIC_NUMBER_MAP.contains_key(input[4*i]))
                .map(|i| {
                    match (input[4*i].parse::<String>(),
                           input[4*i+1].parse::<f64>(),
                           input[4*i+2].parse::<f64>(),
                           input[4*i+3].trim_end_matches(";").parse::<f64>()) {
                        (Ok(a), Ok(b), Ok(c), Ok(d)) => Ok((a,b,c,d)),
                        _ => Err("It appears your string is no good :("),
                    }})
                    .collect::<Result<Vec<(String,f64,f64,f64)>,&str>>()?;
        Ok(())
    }

    fn gaussian_int(&self, n: f64, alpha: f64) -> f64 {
        let n1 = (n+1_f64) / 2_f64;
        gamma(n1) / (2_f64 * alpha.powf(n1))
    }
    fn gto_norm(&self, ang_mom: i32, expnt: f64) -> f64 {
        1_f64 / self.gaussian_int((ang_mom as f64)*2_f64+2_f64, 2_f64*expnt).sqrt()
    }

    pub fn build(&mut self) {
        //TODO make this unique atoms (probably with HashSet)
        let atom_strs = self.atoms.iter().map(|x| x.0.clone()).collect::<Vec<String>>();
        let atom_vec = parse_nwchem::get_basis_info(BASIS_MAP.get(&self.basis).unwrap(),atom_strs).unwrap();
        let c_atm: Vec<Vec<i32>> =  atom_vec.iter().enumerate()
                             .map(|(i,atom)| 
                                    vec![*ATOMIC_NUMBER_MAP.get(atom.get_symb()).unwrap(),
                                    (3*i) as i32 + 20,
                                    0_i32,
                                    0_i32,
                                    0_i32,
                                    0_i32]).collect();
        let mut c_bas: Vec<Vec<Vec<i32>>> = atom_vec.iter().enumerate()
                                    .map(|(i,atom)| 
                                        atom.get_gtos().iter()
                                            .map(|gto| 
                                                vec![
                                                    i as i32,
                                                    gto.get_ang_mom(),
                                                    gto.get_num_prim(),
                                                    gto.get_num_cont(),
                                                    0,//zero for now,
                                                    0,
                                                    0,
                                                    0
                                                ]).collect::<Vec<Vec<i32>>>()
                                            ).collect();

        let mut c_env: Vec<f64> = vec![Vec::from(&[0_f64;20]),self.atoms.iter().map(|atom| vec![atom.1,atom.2,atom.3]).flatten().collect()].concat();

        for (p,atom) in atom_vec.iter().enumerate() {
            for (q,gto) in atom.get_gtos().iter().enumerate(){
                //let mut cs = gto.get_orbs().iter().map(|vec| vec[1..].iter().rev().map(|x|*x).collect::<Vec<f64>>()).collect::<Vec<Vec<f64>>>();
                let mut cs = gto.get_orbs().iter().map(|vec| Vec::from(&vec[1..])).collect::<Vec<Vec<f64>>>();
                let mut es = gto.get_orbs().iter().map(|vec| vec[0]).collect::<Vec<f64>>();
                cs = cs.iter().enumerate().map(|(i,vec)| {
                            let temp = self.gto_norm(gto.get_ang_mom(),es[i]);
                            vec.iter().map(|val| val*temp).collect::<Vec<f64>>()
                            }).collect::<Vec<Vec<f64>>>();

                let s1 = (0..cs[0].len()).map(|k| 
                    iproduct!(0..es.len(),0..es.len())
                        .fold(0_f64,|acc,(i,j)| 
                            acc + self.gaussian_int((gto.get_ang_mom()*2+2).into(),&es[i]+&es[j]) 
                            * &cs[j][k] * &cs[i][k]))
                     .collect::<Vec<f64>>();

                let mut norm = cs.iter()
                             .map(|vec| 
                                 vec.iter().enumerate()
                                    .map(|(i,x)| x/s1[i].sqrt()).collect::<Vec<f64>>())
                             .collect::<Vec<Vec<f64>>>();
                c_bas[p][q][5] = (c_env.len()) as i32;
                c_env.extend(es);
                c_bas[p][q][6] = (c_env.len()) as i32;
                let mut norm_flat: Vec<f64> = Vec::new();
                for i in 0..norm[0].len() {
                    for j in 0..norm.len(){
                        norm_flat.push(norm[j][i]);
                    }
                }
                //c_env.extend(norm.into_iter().flat_map(|vec| vec.into_iter()).collect::<Vec<f64>>())
                c_env.extend(norm_flat)

            }
        }
        //dbg!(&c_env);
        //dbg!(&c_env.iter().fold(0_f64,|acc,x| acc+x.abs()));
        let c_bas = c_bas.concat();
        self.c_dat.initial_r2c(&c_atm,c_atm.len() as i32,&c_bas,c_bas.len() as i32, c_env);
    }

    pub fn get_int1e_ovlp(&mut self) -> Array2<f64>{
        let nbas = self.c_dat.get_nbas();
        let mut vec: Vec<usize> = Vec::new();
        let shls = (0..nbas as i32).map(|i| self.c_dat.cint_cgto_rust(i) as usize).collect::<Vec<usize>>();
        vec.push(0);
        for len in &shls {
            vec.push(vec.last().unwrap()+*len);
        }
        let nao = *vec.last().unwrap();
        let mut out = Array::zeros(nao*nao);
        for i in 0..nbas{
            for j in 0..nbas{
                let s1 = self.c_dat.int1e_ovlp(i as i32, j as i32);
                let indices = iproduct!(vec[j]..vec[j+1],vec[i]..vec[i+1]).collect::<Vec<(usize,usize)>>();
                for (index,item) in s1.iter().enumerate() {
                    out[indices[index].0*nao+indices[index].1] = *item
                }
            }
        }
        out.into_shape_with_order(((nao,nao),Order::ColumnMajor)).unwrap()
    }

    pub fn get_int1e_nuc(&mut self) -> Array2<f64>{
        let nbas = self.c_dat.get_nbas();
        let mut vec: Vec<usize> = Vec::new();
        let shls = (0..nbas as i32).map(|i| self.c_dat.cint_cgto_rust(i) as usize).collect::<Vec<usize>>();
        vec.push(0);
        for len in &shls {
            vec.push(vec.last().unwrap()+*len);
        }
        let nao = *vec.last().unwrap();
        let mut out = Array::zeros(nao*nao);
        for i in 0..nbas{
            for j in 0..nbas{
                let s1 = self.c_dat.int1e_nuc(i as i32, j as i32);
                let indices = iproduct!(vec[j]..vec[j+1],vec[i]..vec[i+1]).collect::<Vec<(usize,usize)>>();
                for (index,item) in s1.iter().enumerate() {
                    out[indices[index].0*nao+indices[index].1] = *item
                }
            }
        }
        out.into_shape_with_order(((nao,nao),Order::ColumnMajor)).unwrap()
    }

    pub fn get_int1e_kin(&mut self) -> Array2<f64>{
        let nbas = self.c_dat.get_nbas();
        let mut vec: Vec<usize> = Vec::new();
        let shls = (0..nbas as i32).map(|i| self.c_dat.cint_cgto_rust(i) as usize).collect::<Vec<usize>>();
        vec.push(0);
        for len in &shls {
            vec.push(vec.last().unwrap()+*len);
        }
        let nao = *vec.last().unwrap();
        let mut out = Array::zeros(nao*nao);
        for i in 0..nbas{
            for j in 0..nbas{
                let s1 = self.c_dat.int1e_kin(i as i32, j as i32);
                let indices = iproduct!(vec[j]..vec[j+1],vec[i]..vec[i+1]).collect::<Vec<(usize,usize)>>();
                for (index,item) in s1.iter().enumerate() {
                    out[indices[index].0*nao+indices[index].1] = *item
                }
            }
        }
        out.into_shape_with_order(((nao,nao),Order::ColumnMajor)).unwrap()
    }

    pub fn get_int2e(&mut self) -> Array4<f64>{
        let nbas = self.c_dat.get_nbas();
        let shls = (0..nbas as i32).map(|i| self.c_dat.cint_cgto_rust(i) as usize).collect::<Vec<usize>>();
        let mut vec: Vec<usize> = Vec::new();
        vec.push(0);
        for len in &shls {
            vec.push(vec.last().unwrap()+*len);
        }
        let nao = *vec.last().unwrap();
        let mut out = Array::zeros(nao*nao*nao*nao);
        for i in 0..nbas {
            for j in 0..nbas {
                for k in 0..nbas {
                    for l in 0..nbas {
                        let s1 = self.c_dat.int2e(i as i32, j as i32, k as i32, l as i32); //Array::from_shape_vec((shls[i],shls[k],shls[j],shls[l]),self.c_dat.int2e(i as i32, j as i32, k as i32, l as i32)).unwrap();
                        let indices = iproduct!(vec[l]..vec[l+1],vec[k]..vec[k+1],vec[j]..vec[j+1],vec[i]..vec[i+1]).collect::<Vec<(usize,usize,usize,usize)>>();
                        for (index,item) in s1.iter().enumerate() {
                            out[indices[index].0*nao.pow(3)+indices[index].1*nao.pow(2)+indices[index].2*nao.pow(1)+indices[index].3] = *item;
                            
                        }
                        //println!("cat");
                        //out.slice_mut(s![vec[i]..vec[i+1],vec[k]..vec[k+1],vec[j]..vec[j+1],vec[l]..vec[l+1]]).assign(&s1);
                    }
                }
            }
        }
        out.into_shape_with_order(((nao,nao,nao,nao),Order::ColumnMajor)).unwrap()
    }
}


#[cfg(test)]
mod tests {
    use ndarray_npy::ReadNpyExt;
    use std::fs::File;
    use super::*;
    #[test]
    fn basis_success() {
        let mut mol = Mol::new();
        assert!(mol.set_basis("sto-3g").is_ok());
        assert!(mol.set_basis("sto-3g!^%&(&^%$").is_ok());
        assert!(mol.set_basis("sto3g     ").is_ok());
        assert!(mol.set_basis("     s   to-3g").is_ok());
        assert!(mol.set_basis("     STO3g").is_ok());
    }

    #[test]
    fn basis_fails() {
        let mut mol = Mol::new();
        assert!(mol.set_basis("sto-3g*").is_err());
        assert!(mol.set_basis("sto-3g!^%&*(*&^%$a").is_err());
        assert!(mol.set_basis("sto3g  b   ").is_err());
        assert!(mol.set_basis("  t   sto-3g").is_err());
        assert!(mol.set_basis("     STO3g++").is_err());
    }

    #[test]
    fn atoms_success() {
        let mut mol = Mol::new();
        assert!(
            mol.set_atoms_str(
            "H 0 0.000000000000 .0;
             He .8 1.0 2.0"
            ).is_ok());
        let test_atoms = Vec::from([("H".to_string(), 0_f64, 0_f64, 0_f64),("He".to_string(), 0.8_f64, 1_f64, 2_f64)]);
        for i in 0..2 {
            assert!(test_atoms[i] == mol.atoms[i])
        }
        let mut mol = Mol::new();
        assert!(
            mol.set_atoms_str(
            "H 0 0.000000000000 .0
             He .8 1.0 2.0"
            ).is_ok());
        let test_atoms = Vec::from([("H".to_string(), 0_f64, 0_f64, 0_f64),("He".to_string(), 0.8_f64, 1_f64, 2_f64)]);
        for i in 0..2 {
            assert!(test_atoms[i] == mol.atoms[i])
        }
    }

    #[test]
    fn test_int1e_ovlp() {
        fn run(coords: &str, basis: &str, file: &str) {
            let mut mol = Mol::new();
            let _ = mol.set_atoms_str(coords);
            let _ = mol.set_basis(basis);
            let _ = mol.build();
            let ovl = mol.get_int1e_ovlp();
            let reader = File::open(format!("{}{}",env!("CARGO_MANIFEST_DIR"),file)).unwrap();
            let temp = Array2::<f64>::read_npy(reader).unwrap();
            dbg!(&ovl.shape());
            dbg!(&temp.shape());
            let thing = ovl-temp;
            assert!(thing.abs().sum()<1e-10);
        }
        run("H 0 0 0; H 1 0 0","sto-3g","/test_data/HH-int1e-ovlp.npy");
        run("H 0 0 0; H 1 0 0","cc-pvdz","/test_data/HH-ccpvdz-int1e-ovlp.npy");
        run("H 0 0 0; Li 1 0 0","sto-3g","/test_data/HLi-int1e-ovlp.npy");
        run("H 0 0 0; Li 1 0 0","cc-pvdz","/test_data/HLi-ccpvdz-int1e-ovlp.npy");
        run("Ti 0 0 0","sto-3g","/test_data/Ti-int1e-ovlp.npy");
        run("Ti 0 0 0","cc-pvdz","/test_data/Ti-ccpvdz-int1e-ovlp.npy");
        run("Ge 0 0 0; As 1 0 0; Br .5 .5 0","cc-pvdz","/test_data/GeAsBr-ccpvdz-int1e-ovlp.npy");
    }

    #[test]
    fn test_int1e_nuc() {
        fn run(coords: &str, basis: &str, file: &str) {
            let mut mol = Mol::new();
            let _ = mol.set_atoms_str(coords);
            let _ = mol.set_basis(basis);
            let _ = mol.build();
            let ovl = mol.get_int1e_nuc();
            let reader = File::open(format!("{}{}",env!("CARGO_MANIFEST_DIR"),file)).unwrap();
            let temp = Array2::<f64>::read_npy(reader).unwrap();
            let thing = ovl-temp;
            assert!(thing.abs().sum()<1e-10);
        }
        run("H 0 0 0; H 1 0 0","sto-3g","/test_data/HH-int1e-nuc.npy");
        run("H 0 0 0; H 1 0 0","cc-pvdz","/test_data/HH-ccpvdz-int1e-nuc.npy");
        run("H 0 0 0; Li 1 0 0","sto-3g","/test_data/HLi-int1e-nuc.npy");
        run("H 0 0 0; Li 1 0 0","cc-pvdz","/test_data/HLi-ccpvdz-int1e-nuc.npy");
        run("Ti 0 0 0","sto-3g","/test_data/Ti-int1e-nuc.npy");
        run("Ti 0 0 0","cc-pvdz","/test_data/Ti-ccpvdz-int1e-nuc.npy");
        run("Ge 0 0 0; As 1 0 0; Br .5 .5 0","cc-pvdz","/test_data/GeAsBr-ccpvdz-int1e-nuc.npy");
    }
    #[test]
    fn test_int1e_kin() {
        fn run(coords: &str, basis: &str, file: &str) {
            let mut mol = Mol::new();
            let _ = mol.set_atoms_str(coords);
            let _ = mol.set_basis(basis);
            let _ = mol.build();
            let ovl = mol.get_int1e_kin();
            let reader = File::open(format!("{}{}",env!("CARGO_MANIFEST_DIR"),file)).unwrap();
            let temp = Array2::<f64>::read_npy(reader).unwrap();
            let thing = ovl-temp;
            assert!(thing.abs().sum()<1e-10);
        }
        run("H 0 0 0; H 1 0 0","sto-3g","/test_data/HH-int1e-kin.npy");
        run("H 0 0 0; H 1 0 0","cc-pvdz","/test_data/HH-ccpvdz-int1e-kin.npy");
        run("H 0 0 0; Li 1 0 0","sto-3g","/test_data/HLi-int1e-kin.npy");
        run("H 0 0 0; Li 1 0 0","cc-pvdz","/test_data/HLi-ccpvdz-int1e-kin.npy");
        run("Ti 0 0 0","sto-3g","/test_data/Ti-int1e-kin.npy");
        run("Ti 0 0 0","cc-pvdz","/test_data/Ti-ccpvdz-int1e-kin.npy");
        run("Ge 0 0 0; As 1 0 0; Br .5 .5 0","cc-pvdz","/test_data/GeAsBr-ccpvdz-int1e-kin.npy");
    }

    #[test]
    fn test_int2e() {
        fn run(coords: &str, basis: &str, file: &str) {
            let mut mol = Mol::new();
            let _ = mol.set_atoms_str(coords);
            let _ = mol.set_basis(basis);
            let _ = mol.build();
            let ovl = mol.get_int2e();
            let reader = File::open(format!("{}{}",env!("CARGO_MANIFEST_DIR"),file)).unwrap();
            let temp = Array4::<f64>::read_npy(reader).unwrap();
            let thing = ovl-temp;
            assert!(thing.abs().sum()<1e-9);
        }
        run("H 0 0 0; H 1 0 0","sto-3g","/test_data/HH-int2e.npy");
        run("H 0 0 0; H 1 0 0","cc-pvdz","/test_data/HH-ccpvdz-int2e.npy");
        run("H 0 0 0; Li 1 0 0","sto-3g","/test_data/HLi-int2e.npy");
        run("H 0 0 0; Li 1 0 0","cc-pvdz","/test_data/HLi-ccpvdz-int2e.npy");
        run("Ti 0 0 0","sto-3g","/test_data/Ti-int2e.npy");
        run("Ti 0 0 0","cc-pvdz","/test_data/Ti-ccpvdz-int2e.npy");
        run("Ge 0 0 0; As 1 0 0; Br .5 .5 0","cc-pvdz","/test_data/GeAsBr-ccpvdz-int2e.npy");
    }
}



