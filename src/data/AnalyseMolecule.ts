/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MolUtil.ts'/>
///<reference path='../../../WebMolKit/src/data/Chemistry.ts'/>
///<reference path='../../../WebMolKit/src/data/BondArtifact.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Find things wrong with a molecule, pertaining to having an appropriate valence.
*/

const VALENCES:{[id:string] : number[]} =
{
	'H': [1],
	'C': [4],
	'N': [3],
	'O': [2],
	'Si': [4],
	'P': [3, 5],
	'S': [2, 4, 6],
	'F': [1],
	'Cl': [1, 3, 5, 7],
	'Br': [1, 3, 5, 7],
	'I': [1, 3, 5, 7],
};

const OXSTATES:{[id:string] : number[]} =
{
};

/*
public static ELEMENTS =
[
	null,
	'H',                                                                                 'He',
	'Li','Be',                                                  'B', 'C', 'N', 'O', 'F', 'Ne',
	'Na','Mg',                                                  'Al','Si','P', 'S', 'Cl','Ar',
	'K', 'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
	'Rb','Sr','Y', 'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I', 'Xe',
	'Cs','Ba',
			'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
					'Lu','Hf','Ta','W', 'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
	'Fr','Ra',
			'Ac','Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
					'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn'
];
*/

export const enum AnalyseMoleculeType
{
	BadValence = 'badvalence', // one of the 'boring' atoms has an invalid valence
	OddOxState = 'oxstate', // one of the 'metals' has an oxidation state that's suspicious
}

export interface AnalyseMoleculeResult
{
	type:AnalyseMoleculeType;
	atom?:number;
	bond?:number;
	value?:number;
}

export class AnalyseMolecule
{
	public results:AnalyseMoleculeResult[] = [];

	// ------------ public methods ------------

	// note that the molecule will be modified if necessary
	constructor(public mol:Molecule)
	{
	}

	public perform():AnalyseMoleculeResult[]
	{
		let mol = this.mol;
		MolUtil.expandAbbrevs(mol, true);

		// look for valence/oxidation state oddities
		for (let n = 1; n <= mol.numAtoms; n++)
		{
			let extra = mol.atomExtra(n), extraMod = false;
			for (let i = extra.length - 1; i >= 0; i--) 
				if (extra[i].startsWith(BONDARTIFACT_EXTRA_RESPATH) || extra[i].startsWith(BONDARTIFACT_EXTRA_RESRING) || extra[i].startsWith(BONDARTIFACT_EXTRA_ARENE))
			{
				extra.splice(i, 1);
				extraMod = true;
			}
			if (extraMod) mol.setAtomExtra(n, extra);

			if (mol.atomAdjCount(n) == 0) continue; // no point in making valence judgments for isolated atoms (?)

			let el = mol.atomElement(n);
			let wantVal = VALENCES[el], wantOx = OXSTATES[el];
			if (!wantVal && !wantOx) continue;

			let hyd = mol.atomHydrogens(n);
			let valence = -mol.atomCharge(n) + mol.atomUnpaired(n) + hyd;
			let oxstate = mol.atomCharge(n) + hyd;
			for (let b of mol.atomAdjBonds(n))
			{
				let o = mol.bondOrder(b);
				valence += o;
				oxstate += o % 2;
			}
			if (wantVal && wantVal.indexOf(valence) < 0)
				this.results.push({'type': AnalyseMoleculeType.BadValence, 'atom': n, 'value': valence});
			if (wantOx && wantOx.indexOf(oxstate) < 0)
				this.results.push({'type': AnalyseMoleculeType.OddOxState, 'atom': n, 'value': oxstate});
		}

		return this.results;
	}

	// ------------ private methods ------------


}

/* EOF */ }