use r_libcint as cint;
use croc_parse::{parse_nwchem, BASIS_MAP};
use statrs::function::gamma::gamma;
use itertools::{Itertools,iproduct};
use crate::ATOMIC_NUMBER_MAP;

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
                            * &cs[j][k] * & cs[i][k]))
                     .collect::<Vec<f64>>();

                let mut norm = cs.iter()
                             .map(|vec| 
                                 vec.iter().enumerate()
                                    .map(|(i,x)| x*s1[i]).collect::<Vec<f64>>())
                             .collect::<Vec<Vec<f64>>>();
                c_bas[p][q][5] = (c_env.len()) as i32;
                c_env.extend(es);
                c_bas[p][q][6] = (c_env.len()) as i32;
                c_env.extend(norm.into_iter().flatten().collect::<Vec<f64>>())
            }
        }
        let c_bas = c_bas.concat();
        self.c_dat.initial_r2c(&c_atm,c_atm.len() as i32,&c_bas,c_bas.len() as i32, c_env);
    }
}


#[cfg(test)]
mod tests {
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

    //#[test]
    //fn test_eri() {
    //    let mut mol = Mol::new();
    //    mol.set_atoms_str(
    //        "H 0 0.000000000000 .0
    //         H 1 0 0");
    //    mol.set_basis("sto-3g");
    //    mol.build();
    //    //dbg!(&mol.c_dat);
    //    //for i in 0..2 {
    //    //    for j in 0..2 {
    //    //        for k in 0..2 {
    //    //            for l in 0..2 {
    //    //                println!("{:?}",mol.c_dat.int2e(i,j,k,l)[0]);
    //    //            }
    //    //        }
    //    //    }
    //    //}
    //}
}



