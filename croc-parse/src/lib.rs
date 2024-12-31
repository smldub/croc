pub mod parse_nwchem;
use phf::{phf_map, Map};
pub static BASIS_MAP: Map::<&'static str,&'static str> = phf_map! {
    "sto3g" => "sto-3g.dat",
    "ccpvdz" => "cc-pvdz.dat",
    "ccpvtz" => "cc-pvtz.dat",
    "augccpvdz" => "aug-cc-pvdz.dat",
};



