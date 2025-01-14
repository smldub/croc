use r_libcint as cint;
use croc_parse::{parse_nwchem, BASIS_MAP};
use statrs::function::gamma::gamma;
use itertools::{Itertools,iproduct};
use crate::ATOMIC_NUMBER_MAP;
use ndarray::{array,Array,Array2,Array4, s, Order};

const BASIS_EXCEPTIONS: [char;2] = ['*','+'];

/// Builder of Mol Object that allows for specification of basis and atom coordinates
pub struct MolBuilder {
    /// The jth atom's name and coords is stored as 
    /// atoms[j] = (name,(x,y,z))
    pub atoms: String,//Vec<(String,f64,f64,f64)>,
    /// Name of basis (see TODO for set of available bases)
    pub basis: String,
    // c_dat contains the correctly formatted data for libcint
    //c_dat: cint::CINTR2CDATA,
}

impl MolBuilder {
    /// Generates the default MolBuilder Object
    pub fn new() -> Self {
        Self{
            atoms: String::new(),
            basis: String::new(),
        }
    }
    
    /// Method for setting basis string
    /// Check BASIS_MAP in lib.rs for available bases
    pub fn set_basis(mut self, basis_name: &str) -> MolBuilder {
        self.basis = String::from(basis_name);
        self
    }

    /// Method for setting atoms coordinates in xyz format i.e.
    /// Format = "a_1 x_1 y_1 z_1; a_2 x_2 y_2 z_2;...."
    /// where a_i is atomic symbol for ith atom and 
    /// x_i,y_i,z_i are cartesian coordinates for ith atom (in bohr)
    pub fn set_atoms_str(mut self, atm_coords: &str) -> MolBuilder {
        self.atoms = String::from(atm_coords);
        self
    }


    /// Internal method for ensuring that basis is in correct format
    /// and is contained in BASIS_MAP
    fn format_basis(&self, basis_name: &str) -> Result<String,&str> {
        let temp = basis_name.trim()
                             .to_lowercase()
                             .chars()
                             .filter(|x| x.is_alphanumeric() | &BASIS_EXCEPTIONS.contains(x))
                             .collect::<String>();
        if BASIS_MAP.contains_key(temp.as_str()) {
            Ok(temp)
        }
        else {
            Err("Basis is not included in std library")
        }
    }

    /// Internal method for ensuring that coordinates have been entered correctly
    /// i.e. xyz format
    fn format_atoms_str(&self, atm_coords: &str) 
        -> Result<Vec<(String,f64,f64,f64)>,&str> {
        let input = atm_coords.trim().split_whitespace().collect::<Vec<&str>>();

        (0..input.len()/4)
                .filter(|i| ATOMIC_NUMBER_MAP.contains_key(input[4*i]))
                .map(|i| {
                    match (input[4*i].parse::<String>(),
                           input[4*i+1].parse::<f64>(),
                           input[4*i+2].parse::<f64>(),
                           input[4*i+3].trim_end_matches(";").parse::<f64>()) {
                        (Ok(a), Ok(b), Ok(c), Ok(d)) => Ok((a,b,c,d)),
                        _ => Err("It appears your atm_coord is no good :("),
                    }})
                    .collect::<Result<Vec<(String,f64,f64,f64)>,&str>>()
    }

    /// TODO
    fn gaussian_int(&self, n: f64, alpha: f64) -> f64 {
        let n1 = (n+1_f64) / 2_f64;
        gamma(n1) / (2_f64 * alpha.powf(n1))
    }

    /// TODO
    fn gto_norm(&self, ang_mom: i32, expnt: f64) -> f64 {
        1_f64 / self.gaussian_int((ang_mom as f64)*2_f64+2_f64, 2_f64*expnt).sqrt()
    }

    /// Build method for Mol object
    pub fn build(&self) -> Result<Mol,&str> {
        // Check that basis and atom info were formatted correctly
        let basis = self.format_basis(self.basis.as_str())?;
        let atoms = self.format_atoms_str(self.atoms.as_str())?;

        //TODO make this unique atoms (probably with HashSet)
        // Parse coordinate information for the atomic symbols
        let atom_strs = atoms.iter().map(|x| x.0.clone()).collect::<Vec<String>>();
        // Get the basis information for the atomic symbols
        let atom_vec = parse_nwchem::get_basis_info(BASIS_MAP.get(&basis).unwrap(),atom_strs).unwrap();
        
        // Create the c_atm array in the format of libcint (see r-libcint for format)
        let c_atm: Vec<Vec<i32>> =  atom_vec.iter().enumerate()
                             .map(|(i,atom)| 
                                    vec![*ATOMIC_NUMBER_MAP.get(atom.get_symb()).unwrap(),
                                    (3*i) as i32 + 20,
                                    0_i32,
                                    0_i32,
                                    0_i32,
                                    0_i32]).collect();

        // Initializing the c_bas array in the format of libcint (see r-libcint for format)
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

        // Initializing the c_env array in the format of libcint (see r-libcint for format)
        let mut c_env: Vec<f64> = vec![Vec::from(&[0_f64;20]),atoms.iter().map(|atom| vec![atom.1,atom.2,atom.3]).flatten().collect()].concat();

        // In this section we normalize our GTOs, put that 
        // information into the c_env array. We also store index
        // information in the c_bas array
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
                c_env.extend(norm_flat)
            }
        }
        // Finally, we can format everything for the creation of the CINTR2CDATA object (I didn't name it)
        let c_bas = c_bas.concat();
        let mut c_dat = cint::CINTR2CDATA::new();
        c_dat.initial_r2c(&c_atm,c_atm.len() as i32,&c_bas,c_bas.len() as i32, c_env);

        // Actually finally, we generate our Mol object (see Mol Struct for description of these
        // objects)
        let nbas = c_dat.get_nbas();
        let shls = (0..nbas as i32).map(|i| c_dat.cint_cgto_rust(i) as usize).collect::<Vec<usize>>();
        let mut dims: Vec<usize> = Vec::with_capacity(shls.len()+1);
        dims.push(0);
        dims.extend(shls.iter().scan(0,|acc, len| {*acc += *len; Some(*acc)} ).collect::<Vec<usize>>());
        let nao = *dims.last().unwrap();
        Ok(Mol {
            atoms: atoms,
            basis: basis,
            c_dat: c_dat,
            nbas: nbas,
            nao: nao,
            shls: shls,
            dims: dims,
        })
    }
}

pub struct Mol {
    /// The jth atom's name and coords is stored as 
    /// atoms[j] = (name,(x,y,z))
    atoms: Vec<(String,f64,f64,f64)>,
    /// Name of basis (see TODO for set of available bases)
    basis: String,
    /// c_dat contains the correctly formatted data for libcint
    c_dat: cint::CINTR2CDATA,
    /// Number of Basis Functions Shells
    nbas: usize,
    /// Number of Atomic Orbitals
    nao: usize,
    /// Vector of shell sizes (i.e. number of atomic orbitals)
    shls: Vec<usize>,
    /// Vector of shell indices (i.e. number of atomic orbitals before a given shell)
    dims: Vec<usize>,
    
}

/// The Mol object is intended to be the immutable generator of electron integrals
impl Mol {
    /// Retrieves the number of atomic orbitals
    pub fn get_nao(&self) -> usize {
        self.nao.clone()
    }

    /// Retrieves the number of basis shells
    pub fn get_nbas(&self) -> usize {
        self.nbas.clone()
    }

    /// Retrieves the atoms and their coordinates (in bohr)
    pub fn get_atoms(&self) -> Vec<(String,f64,f64,f64)> {
        self.atoms.clone()
    }

    /// Generates the overlap integrals 
    /// TODO decide on mathematical representation
    pub fn get_int1e_ovlp(&self) -> Array2<f64>{
        // Generate array for 1e integrals 
        let mut out = Array::zeros((self.nao,self.nao));
        // Iterate through the shells
        for j in 0..self.nbas {
            for i in j..self.nbas {
                // Temporarily store the COLUMNMAJOR ordered integrals between the ith and jth
                // shells
                let s1 = self.c_dat.int1e_ovlp(i as i32, j as i32);
                // Iterate through the indices of the shells, taking advantage of symmetry
                // in ColumnMajor Order (j changes slowest then i)
                for (index, (k,l)) in iproduct!(self.dims[j]..self.dims[j+1],self.dims[i]..self.dims[i+1]).enumerate() {
                    out[[k,l]] = s1[index];
                    out[[l,k]] = s1[index];
                }
            }
        }
        out
    }

    //pub fn get_int1e_ovlp_b(&self) -> Array2<f64>{
    //    // Generate array for 1e integrals 
    //    let mut out = Array::zeros((self.nao,self.nao));
    //    // Iterate through the shells
    //    for j in 0..self.nbas{
    //        for i in j..self.nbas{
    //            // Temporarily store the COLUMNMAJOR ordered integrals between the ith and jth
    //            // shells
    //            let s1 = Array::from_vec(self.c_dat.int1e_ovlp(i as i32, j as i32))
    //                .into_shape_with_order(((self.shls[i],self.shls[j]),Order::ColumnMajor)).unwrap();
    //            // Assign slices (taking advantage of symmetry)
    //            // benchmarking vs indexing (see method below) shows that calls to underlying 
    //            // C-program are pretty inconsistent in timing making it hard to benchmark
    //            out.slice_mut(s![self.dims[i]..self.dims[i+1],self.dims[j]..self.dims[j+1]]).assign(&s1);
    //            out.slice_mut(s![self.dims[j]..self.dims[j+1],self.dims[i]..self.dims[i+1]]).assign(&s1.t());
    //            }
    //        }
    //    }
    //    out
    //}

    /// Generates the nuclear 1-e integrals
    /// TODO decide on mathematical representation
    pub fn get_int1e_nuc(&self) -> Array2<f64>{
        // Generate array for 1e integrals 
        let mut out = Array::zeros((self.nao,self.nao));
        // Iterate through the shells
        for j in 0..self.nbas {
            for i in j..self.nbas {
                // Temporarily store the COLUMNMAJOR ordered integrals between the ith and jth
                // shells
                let s1 = self.c_dat.int1e_nuc(i as i32, j as i32);
                // Iterate through the indices of the shells, taking advantage of symmetry
                // in ColumnMajor Order (j changes slowest then i)
                for (index, (k,l)) in iproduct!(self.dims[j]..self.dims[j+1],self.dims[i]..self.dims[i+1]).enumerate() {
                    out[[k,l]] = s1[index];
                    out[[l,k]] = s1[index];
                }
            }
        }
        out
    }

    /// Generates the kinetic 1-e integrals
    /// TODO decide on mathematical representation
    pub fn get_int1e_kin(&self) -> Array2<f64>{
        // Generate array for 1e integrals 
        let mut out = Array::zeros((self.nao,self.nao));
        // Iterate through the shells
        for j in 0..self.nbas {
            for i in j..self.nbas {
                // Temporarily store the COLUMNMAJOR ordered integrals between the ith and jth
                // shells
                let s1 = self.c_dat.int1e_kin(i as i32, j as i32);
                // Iterate through the indices of the shells, taking advantage of symmetry
                // in ColumnMajor Order (j changes slowest then i)
                for (index, (k,l)) in iproduct!(self.dims[j]..self.dims[j+1],self.dims[i]..self.dims[i+1]).enumerate() {
                    out[[k,l]] = s1[index];
                    out[[l,k]] = s1[index];
                }
            }
        }
        out
    }

    /// Generates the Electron Repulsion Integrals (ERIs) with indices matching chemists notation
    /// (ij|kl)
    pub fn get_int2e(&self) -> Array4<f64>{
        // Generate array for 2e integrals in chemistry notation (ij|kl)
        let mut out = Array::zeros((self.nao,self.nao,self.nao,self.nao));
        // We can take advantage of 8 fold symmetry in ERIs
        // ijkl = jikl, therefore we only need combinations (with replacement) of ij
        // the same applies to kl. 
        // Additionally we know that ijkl = klij therefore we only need combinations (with
        // replacement) of pairs of ij and kl
        let comb = (0..self.nbas).combinations_with_replacement(2).combinations_with_replacement(2);
        for indices in comb {
            // Writing out indices for easy of use (should not effect speed)
            let (i,j,k,l) = (indices[0][0],indices[0][1],indices[1][0],indices[1][1]);
            // Getting eris for shells i,j,k,l
            let s1 = self.c_dat.int2e(i as i32, j as i32, k as i32, l as i32);
            // Writing out iterator over indices in the shells i,j,k,l
            // in ColumnMajor Order (l changes slowest then k then j ...)
            let iter = iproduct!(self.dims[l]..self.dims[l+1],
                                    self.dims[k]..self.dims[k+1],
                                    self.dims[j]..self.dims[j+1],
                                    self.dims[i]..self.dims[i+1]);
            // Finally, assigning values to eri array taking advantage of 8-fold symmetry
            for (index,(i1,i2,i3,i4)) in iter.enumerate() {
                out[[i1,i2,i3,i4]] = s1[index];
                out[[i2,i1,i3,i4]] = s1[index];
                out[[i1,i2,i4,i3]] = s1[index];
                out[[i2,i1,i4,i3]] = s1[index];
                out[[i3,i4,i1,i2]] = s1[index];
                out[[i4,i3,i1,i2]] = s1[index];
                out[[i3,i4,i2,i1]] = s1[index];
                out[[i4,i3,i2,i1]] = s1[index];
            }
        }
        out
    }
}


#[cfg(test)]
mod tests {
    use ndarray_npy::ReadNpyExt;
    use std::fs::File;
    use super::*;
    #[test]
    fn basis_success() {
        let mut mol = MolBuilder::new();
        assert!(mol.format_basis("sto-3g").is_ok());
        assert!(mol.format_basis("sto-3g!^%&(&^%$").is_ok());
        assert!(mol.format_basis("sto3g     ").is_ok());
        assert!(mol.format_basis("     s   to-3g").is_ok());
        assert!(mol.format_basis("     STO3g").is_ok());
    }

    #[test]
    fn basis_fails() {
        let mut mol = MolBuilder::new();
        assert!(mol.format_basis("sto-3g*").is_err());
        assert!(mol.format_basis("sto-3g!^%&*(*&^%$a").is_err());
        assert!(mol.format_basis("sto3g  b   ").is_err());
        assert!(mol.format_basis("  t   sto-3g").is_err());
        assert!(mol.format_basis("     STO3g++").is_err());
    }

    #[test]
    fn atoms_success() {
        let mut mol = MolBuilder::new();
        let result_atoms = mol.format_atoms_str(
            "H 0 0.000000000000 .0;
             He .8 1.0 2.0"
            );
        assert!(result_atoms.is_ok());
        let result_atoms = result_atoms.unwrap();
        let test_atoms = Vec::from([("H".to_string(), 0_f64, 0_f64, 0_f64),("He".to_string(), 0.8_f64, 1_f64, 2_f64)]);
        for i in 0..2 {
            assert!(test_atoms[i] == result_atoms[i])
        }
        let mut mol = MolBuilder::new();
        let result_atoms = mol.format_atoms_str(
            "H 0 0.000000000000 .0;
             He .8 1.0 2.0"
            );
        assert!(result_atoms.is_ok());
        let result_atoms = result_atoms.unwrap();
        //let test_atoms = Vec::from([("H".to_string(), 0_f64, 0_f64, 0_f64),("He".to_string(), 0.8_f64, 1_f64, 2_f64)]);
        for i in 0..2 {
            assert!(test_atoms[i] == result_atoms[i])
        }
    }

    #[test]
    fn test_int1e_ovlp() {
        fn run(coords: &str, basis: &str, file: &str) {
            let mut mol = MolBuilder::new()
                            .set_atoms_str(coords)
                            .set_basis(basis)
                            .build().unwrap();
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
        //assert!(1==2);
        run("H 0 0 0; Li 1 0 0","cc-pvdz","/test_data/HLi-ccpvdz-int1e-ovlp.npy");
        run("Ti 0 0 0","sto-3g","/test_data/Ti-int1e-ovlp.npy");
        run("Ti 0 0 0","cc-pvdz","/test_data/Ti-ccpvdz-int1e-ovlp.npy");
        run("Ge 0 0 0; As 1 0 0; Br .5 .5 0","cc-pvdz","/test_data/GeAsBr-ccpvdz-int1e-ovlp.npy");

    }

    #[test]
    fn test_int1e_nuc() {
        fn run(coords: &str, basis: &str, file: &str) {
            let mut mol = MolBuilder::new()
                          .set_atoms_str(coords)
                          .set_basis(basis)
                          .build().unwrap();
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
            let mut mol = MolBuilder::new()
                          .set_atoms_str(coords)
                          .set_basis(basis)
                          .build().unwrap();
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
            let mut mol = MolBuilder::new()
                          .set_atoms_str(coords)
                          .set_basis(basis)
                          .build().unwrap();
            let ovl = mol.get_int2e();
            let reader = File::open(format!("{}{}",env!("CARGO_MANIFEST_DIR"),file)).unwrap();
            let temp = Array4::<f64>::read_npy(reader).unwrap();
            dbg!(&ovl);
            dbg!(&temp);
            let thing = ovl-temp;
            assert!(thing.abs().sum()<1e-9);
        }
        run("H 0 0 0; H 1 0 0","sto-3g","/test_data/HH-int2e.npy");
        run("H 0 0 0; H 1 0 0","cc-pvdz","/test_data/HH-ccpvdz-int2e.npy");
        run("H 0 0 0; Li 1 0 0","sto-3g","/test_data/HLi-int2e.npy");
        run("H 0 0 0; Li 1 0 0","cc-pvdz","/test_data/HLi-ccpvdz-int2e.npy");
        run("Ti 0 0 0","sto-3g","/test_data/Ti-int2e.npy");
        run("Ti 0 0 0","cc-pvdz","/test_data/Ti-ccpvdz-int2e.npy");
    }
}



