#![allow(non_upper_case_globals)]
#![allow(unused)]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
use std::os::raw::c_int;
use std::mem::ManuallyDrop;

#[derive(Debug)]
pub enum CintType {
   Spheric,
   Cartesian,
   //Spinor,  // Not yet included
}
#[derive(Debug)]
pub struct CINTR2CDATA {
    c_atm: (*mut i32, usize, usize),
    c_dims: (*mut i32, usize, usize),
    c_bas: (*mut i32, usize, usize),
    c_env: (*mut f64, usize, usize),
    c_cache:(*mut f64, usize, usize),
    c_nbas: c_int,
    c_natm: c_int,
    c_opt: (*mut CINTOpt, usize, usize),
    cint_type: CintType,
}

impl CINTR2CDATA {
    /// create a new, empty CINTR2CDATA.
    /// //TODO: this is a mishmash of description of inputs for initial_r2c and actual libcint
    /// format
    ///
    /// Basis is stored in libcint format
    /// For jth basis, the eight slots are
    /// c_bas[j][0] = 0-based index of corresponding atom
    /// c_bas[j][1] = angular momentum
    /// c_bas[j][2] = number of primitive GTO
    /// c_bas[j][3] = number of contracted GTO
    /// c_bas[j][4] = kappa for spinor GTO TODO: include more info and understand better
    /// c_bas[j][5] = env offset to save exponents of primitive GTOS (see env definition for
    /// location)
    /// c_bas[j][6] = env offset to save column-major contraction coefficients (see env definition
    /// for location)
    /// c_bas[j][7] = unused
    ///
    /// Env array stores information about coordinates of molecules as as well as 
    /// GTO coefficients (First the primative values, then the normalized contracted
    /// coeffecients)
    /// HOWEVER for some unknown reason the first 20 slots in the env array MUST be 0_f64
    /// c_env[0..20] = 0
    /// c_env[20 + 3*i]: ith atom x coord (in bohr)
    /// c_env[20 + 3*i+1]: ith atom y coord (in bohr) 
    /// c_env[20 + 3*i+2]: ith atom z coord (in bohr)
    /// Let x = 20+3*(atomNum+1) or in otherwords, we will be describing
    /// the data stored after the coordinates
    /// c_env[x .. x + numPrimGTO]: The primative GTO coeff
    /// c_env[x+numPrimGTO..x+numPrimGto+numContGTO]: Contracted GTO coeff
    /// This pattern will repeat for each GTO
    /// because the number of prim and cont coefficients will change
    /// by GTO it is necessary to store "env offset" (index of start of cont and prim
    /// coeffs) in basis array for each GTO
    ///
    /// c_atm[i][0]: ith atom atomic number
    /// c_atm[i][1]: ith atom coord offset (normally 20+i*3, see env definition)
    /// c_atm[i][2]: TODO
    /// c_atm[i][3]: TODO
    /// c_atm[i][4]: TODO
    /// c_atm[i][5]: TODO

    pub fn new() -> CINTR2CDATA {
        CINTR2CDATA { 
            // name: (ptr, len, capacity)
            c_atm:  (std::ptr::null_mut::<i32>(), 0, 0),
            c_dims: (std::ptr::null_mut::<i32>(), 0, 0),
            c_bas:  (std::ptr::null_mut::<i32>(), 0, 0),
            c_env:  (std::ptr::null_mut::<f64>(), 0, 0),
            c_opt:  (std::ptr::null_mut::<CINTOpt>(), 0,0),
            c_cache:(std::ptr::null_mut::<f64>(), 0, 0),
            c_nbas: 0 as c_int,
            c_natm: 0 as c_int,
            cint_type: CintType::Cartesian,
            }
    }
        pub fn set_cint_type(&mut self, ctype: CintType) {
        self.cint_type = ctype;
    }

    pub fn initial_r2c(&mut self, 
                    atm: &Vec<Vec<i32>>, natm:i32, 
                    bas: &Vec<Vec<i32>>, nbas:i32, 
                    mut env: Vec<f64>) {
        env.shrink_to_fit();
        let mut env = ManuallyDrop::new(env);
        self.c_env = (env.as_mut_ptr(), env.len(), env.capacity());

        let mut bas_f: Vec<i32> = bas.concat();
        bas_f.shrink_to_fit();
        let mut bas_f = ManuallyDrop::new(bas_f);
        self.c_bas = (bas_f.as_mut_ptr(), bas_f.len(), bas_f.capacity());

        let mut atm_f:  Vec<i32> = atm.concat();
        atm_f.shrink_to_fit();
        let mut atm_f = ManuallyDrop::new(atm_f);
        self.c_atm = (atm_f.as_mut_ptr(), atm_f.len(), atm_f.capacity());

        self.c_natm = natm as c_int;
        self.c_nbas = nbas as c_int;
    }
    // Are these functions still necessary?
    //pub fn final_c2r(&mut self) {
    //    println!("Transfer the ownership of the raw pointers in CINTR2CDATA to Rust");
    //    unsafe {
    //        let r_atm = Vec::from_raw_parts(self.c_atm.0,self.c_atm.1,self.c_atm.2);
    //        let r_bas = Vec::from_raw_parts(self.c_bas.0,self.c_bas.1,self.c_bas.2);
    //        let r_env = Vec::from_raw_parts(self.c_env.0,self.c_env.1,self.c_env.2);
    //    }
    //    self.cint_del_optimizer_rust();
    //}
    //pub fn cint_del_optimizer_rust(&mut self) {
    //    unsafe{
    //        CINTdel_optimizer(&mut self.c_opt.0);
    //    }
    //}
    pub fn cint2e_optimizer_rust(&mut self){
        //self.cint_del_optimizer_rust();
        //self.cint_init_2e_optimizer_rust();
        unsafe {
            int2e_optimizer(&mut self.c_opt.0, 
                                       self.c_atm.0, self.c_natm, 
                                       self.c_bas.0, self.c_nbas, 
                                       self.c_env.0);
        }
    }
    pub fn int1e_ovlp_optimizer_rust(&mut self){
        //self.cint_del_optimizer_rust();
        //self.cint_init_optimizer_rust();
        unsafe {
            int1e_ovlp_optimizer(&mut self.c_opt.0,
                                       self.c_atm.0, self.c_natm,
                                       self.c_bas.0, self.c_nbas,
                                       self.c_env.0);
        }
    }
    pub fn int1e_nuc_optimizer_rust(&mut self){
        //self.cint_del_optimizer_rust();
        //self.cint_init_optimizer_rust();
        unsafe {
            int1e_nuc_optimizer(&mut self.c_opt.0,
                                       self.c_atm.0, self.c_natm,
                                       self.c_bas.0, self.c_nbas,
                                       self.c_env.0);
        }
    }
    pub fn int1e_kin_optimizer_rust(&mut self){
        //self.cint_del_optimizer_rust();
        //self.cint_init_optimizer_rust();
        unsafe {
            int1e_kin_optimizer(&mut self.c_opt.0,
                                       self.c_atm.0, self.c_natm,
                                       self.c_bas.0, self.c_nbas,
                                       self.c_env.0);
        }
    }
    pub fn cint_cgto_rust(&self, index: i32) -> i32 {
        let mut dim: i32;
        unsafe {
            dim = match self.cint_type {
                CintType::Spheric  => CINTcgto_spheric(index as c_int, self.c_bas.0) as i32,
                CintType::Cartesian=> CINTcgto_cart(index as c_int, self.c_bas.0) as i32,
            };
        }
        dim
    }
    pub fn int2e(&mut self, i:i32,j:i32,k:i32,l:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut dk = self.cint_cgto_rust(k);
        let mut dl = self.cint_cgto_rust(l);
    
        let mut shls: [c_int; 4] = [i as c_int,j as c_int, k as c_int, l as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let c_shls = shls.as_mut_ptr();
    
        let buf_len = (di*dj*dk*dl) as usize;
        let mut buf: Vec<f64> = Vec::with_capacity(buf_len);
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.capacity());

        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => int2e_sph(c_buf,self.c_dims.0, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0,self.c_cache.0,),
                CintType::Cartesian => int2e_cart(c_buf, self.c_dims.0, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0, self.c_cache.0,),
            };
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }
    pub fn int1e_ovlp(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: [c_int; 2] = [i as c_int,j as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let c_shls = shls.as_mut_ptr();
    
        let buf_len = (di*dj) as usize;
        let mut buf: Vec<f64> = Vec::with_capacity(buf_len);
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => int1e_ovlp_sph(
                           c_buf, self.c_dims.0, c_shls,
                             self.c_atm.0, self.c_natm,
                             self.c_bas.0,self.c_nbas,
                             self.c_env.0,
                             self.c_opt.0, self.c_cache.0),
                CintType::Cartesian => int1e_ovlp_cart(
                           c_buf, self.c_dims.0, c_shls,
                             self.c_atm.0, self.c_natm,
                             self.c_bas.0,self.c_nbas,
                             self.c_env.0,
                             self.c_opt.0, self.c_cache.0),
            };
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }
    pub fn int1e_nuc(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: [c_int; 2] = [i as c_int,j as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let c_shls = shls.as_mut_ptr();
    
        let buf_len = (di*dj) as usize;
        let mut buf: Vec<f64> = Vec::with_capacity(buf_len);
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric =>  int1e_nuc_sph(
                           c_buf, self.c_dims.0, c_shls,
                             self.c_atm.0, self.c_natm,
                             self.c_bas.0,self.c_nbas,
                             self.c_env.0,
                             self.c_opt.0, self.c_cache.0,),
                CintType::Cartesian =>  int1e_nuc_cart(
                           c_buf, self.c_dims.0, c_shls,
                             self.c_atm.0, self.c_natm,
                             self.c_bas.0,self.c_nbas,
                             self.c_env.0,
                             self.c_opt.0, self.c_cache.0,),
            };
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }

    pub fn int1e_kin(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);

        let mut shls: [c_int; 2] = [i as c_int,j as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let c_shls = shls.as_mut_ptr();
    
        let buf_len = (di*dj) as usize;
        let mut buf: Vec<f64> = Vec::with_capacity(buf_len);
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.capacity());
   
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => int1e_kin_sph(
                              c_buf, self.c_dims.0, c_shls,
                                self.c_atm.0, self.c_natm,
                                self.c_bas.0,self.c_nbas,
                                self.c_env.0,
                                self.c_opt.0, self.c_cache.0),
                CintType::Cartesian => int1e_kin_cart(
                              c_buf, self.c_dims.0, c_shls,
                                self.c_atm.0, self.c_natm,
                                self.c_bas.0,self.c_nbas,
                                self.c_env.0,
                                self.c_opt.0, self.c_cache.0),
            };
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }    
}
#[cfg(test)]
mod tests {
    use super::*;
    
    pub fn test_init() -> CINTR2CDATA {
        let mut atm: Vec<Vec<i32>> = vec![];
        atm.push(vec![1,0,0,0,0,0]);
        atm.push(vec![4,3,0,0,0,0]);
        let natm = atm.len() as i32;
        let mut bas: Vec<Vec<i32>> = vec![];
        bas.push(vec![0,1,1,1,1,6,7,0]);
        bas.push(vec![0,2,1,1,1,8,9,0]);
        let nbas = bas.len() as i32;
        let mut env: Vec<f64> = vec![0.0,0.0,0.0,0.7,0.0,0.0,1.0,1.0,0.5,1.0];
        let mut cint_data = CINTR2CDATA::new();
        cint_data.initial_r2c(&atm,natm,&bas,nbas,env);
        cint_data
    }

    #[test]
    fn int2e_sph() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Spheric);
        let buf = cint_data.int2e(0,0,0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((0.0588972110870412-v1).abs() < 1.0e-6);
    }

    #[test]
    fn int2e_cart() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Cartesian);
        let buf = cint_data.int2e(0,0,0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((0.0588972110870412-v1).abs() < 1.0e-6);
    }

    #[test]
    fn int1e_ovlp_sph() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Spheric);
        let buf = cint_data.int1e_ovlp(0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((0.35249460111998443-v1).abs() < 1.0e-6);
    }

    #[test]
    fn int1e_ovlp_cart() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Cartesian);
        let buf = cint_data.int1e_ovlp(0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((0.35249460111998443-v1).abs() < 1.0e-6);
    }

    #[test]
    fn int1e_nuc_sph() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Spheric);
        let buf = cint_data.int1e_nuc(0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((1.7824425521540834-v1).abs() < 1.0e-6);
    }

    #[test]
    fn int1e_nuc_cart() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Cartesian);
        let buf = cint_data.int1e_nuc(0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((1.7824425521540834-v1).abs() < 1.0e-6);
    }

    #[test]
    fn int1e_kin_sph() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Spheric);
        let buf = cint_data.int1e_kin(0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((0.8812365027999611-v1).abs() < 1.0e-6);
    }

    #[test]
    fn int1e_kin_cart() {
        let mut cint_data = test_init();
        cint_data.set_cint_type(CintType::Cartesian);
        let buf = cint_data.int1e_kin(0,0);
        let mut v1:f64=0.0;
        let _ = &buf.iter().for_each(|i| {v1 += i.abs()});
        assert!((0.8812365027999611-v1).abs() < 1.0e-6);
    }
}
