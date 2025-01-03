use std::path::Path;
use std::fs;
use nom::bytes::complete::{tag,take_until};
use nom::multi::many0;
use nom::sequence::{separated_pair, tuple, preceded,terminated};
use nom::character::complete::{alpha0, space0, space1, newline};
use nom::number::complete::double;
use nom::IResult;


#[derive(Debug,PartialEq)]
pub struct Atom {
    //Atomic symbol (e.g. He for Helium)
    symb: String,
    //The associated basis functions grouped by angular momentum
    gtos: Vec<GTO>,
}

impl Atom {
    pub fn get_symb(&self)-> &str{
        self.symb.as_str()
    }
    pub fn get_gtos(&self)-> &Vec<GTO>{
        &self.gtos
    }
}

//TODO: potentially derive getters to cut down on boilerplate
#[derive(Debug,PartialEq)]
pub struct GTO {
    //Angular momentum
    ang_mom: i32,
    //Number of primative orbitals (i.e. number of rows)
    num_prim: i32,
    //Number of contracted orbitals (i.e. number of columns-1)
    num_cont: i32,
    //Orbital information Vec::from([[f64;num_cont+1]; num_prim])
    orbs: Vec<Vec<f64>>,
}

const ANG_ARRAY: [&str;7] = ["S","P","D","F","G","H","I"];

impl GTO {
    pub fn get_ang_mom(&self)-> i32{
        self.ang_mom.clone()
    }
    pub fn get_num_prim(&self)-> i32{
        self.num_prim.clone()
    }
    pub fn get_num_cont(&self)-> i32{
        self.num_cont.clone()
    }

    pub fn get_orbs(&self)-> Vec<Vec<f64>>{
        self.orbs.clone()
    }

    pub fn from_parsed_output<'a>(input: ((&'a str, &'a str),Vec<Vec<f64>>)) -> Option<Vec<Self>>{
       if input.1.len()==0{
           return None;
       }

       if (input.0).1=="SP"{
            //some assumptions are that num_cont is always 1
            return Some(Vec::from([
                GTO { 
                    ang_mom: 0,
                    num_prim: input.1.len() as i32,
                    num_cont: 1,
                    orbs: input.1.iter().map(|vec| Vec::from([vec[0],vec[1]])).collect(),},
                GTO {
                    ang_mom: 1,
                    num_prim: input.1.len() as i32,
                    num_cont: 1, 
                    orbs: input.1.iter().map(|vec| Vec::from([vec[0],vec[2]])).collect(),},
            ]));
       }

       if let Some(ang_mom) = &ANG_ARRAY.iter().position(|x| *x==(input.0).1){
           //TODO check to make sure -1 doesn't overflow
           return Some(vec![GTO {
                                ang_mom: *ang_mom as i32,
                                num_prim: input.1.len() as i32, 
                                num_cont: (input.1[0].len()-1) as i32,
                                orbs: input.1}]);
       }
       None
    }
}

fn parse_file_for_atom<'a,'b>(input: &'a str,atom: &'b str) -> IResult<&'a str,Vec<((&'a str, &'a str), Vec<Vec<f64>>)>> {
    //Discard everything before the PRINT statement
    preceded(take_until("#BA"),
    many0(
        preceded(take_until(atom), 
            tuple((
                terminated(separated_pair(tag(atom),space1,alpha0),newline),
                many0(terminated(many0(preceded(space0, double)),tuple((space0,newline))))
            ))
        )
    ))(input)
}

pub fn parse(path: &Path, atom_list: Vec<String>) -> Result<Vec<Atom>, ()> {
    //TODO: how to deal with this error
    let file = fs::read_to_string(path).unwrap();
    let mut atoms: Vec<Atom> = Vec::new();
    for atom in atom_list.iter() {
        match parse_file_for_atom(file.as_str(), atom.as_str()) {
            Ok(a) => {
                let mut t_atom = Atom {
                    symb: String::from(atom),
                    gtos: a.1.into_iter()
                            .map(|x| GTO::from_parsed_output(x).unwrap()).flatten()
                            .collect::<Vec<GTO>>()
                };
                t_atom.gtos.sort_by_key(|gto| gto.ang_mom);
                atoms.push(t_atom);
            },
            _ => println!("Error"),
        }
    }
    Ok(atoms)
}

pub fn get_basis_info(b_name: &str, atom_list: Vec<String>) -> Result<Vec<Atom>,()> {
        let path_string = format!("{}{}{}",env!("CARGO_MANIFEST_DIR"), "/basis_files/",b_name);
        let path = Path::new(path_string.as_str());
        parse(path, atom_list)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_test() {
        let path = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/basis_files/sto-3g.dat"));
        let atom_list = Vec::from([String::from("H")]);
        let atoms = vec![Atom{symb:String::from("H"), gtos: vec![GTO{ang_mom:0,num_prim:3,num_cont:1,orbs:vec![vec![3.42525091,0.15432897],vec![0.62391373,0.53532814],vec![0.16885540,0.44463454]]}]}];
        assert_eq!(atoms,parse(path,atom_list).unwrap());
    }
}      
