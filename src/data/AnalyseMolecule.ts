/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
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

// normal allowed-valence list
const VALENCES:Record<string, number[]> =
{
	'H': [1],
	'C': [4],
	'N': [3],
	'O': [2],
	'Si': [4],
	'P': [3, 5],
	'S': [2, 4, 6],
	'F': [1],
	'Cl': [1],
	'Br': [1],
	'I': [1],
	'B': [3, 5],
	'Ne': [0],
	'Ar': [0],
	'Ge': [4],
	'As': [3],
	'Se': [2, 4, 6],
	'Kr': [0],
	'Te': [2, 4, 6],
	'Xe': [0],
	'Rn': [0],
};

// oxidised valence options: these are only allowed when an atom has oxygens or fluorines
const OXVALENCES:Record<string, number[]> =
{
	'P': [7],
	'Cl': [3, 5, 7],
	'Br': [3, 5, 7],
	'I': [3, 5, 7],
	'As': [5, 7],
	'Xe': [2, 4, 6],
};

// valence options available only in cationic form: lower than usual valences are seen in some cases (e.g. carbocations)
const CATVALENCES:Record<string, number[]> =
{
	'C': [2],
	'Si': [2],
	'Al': [1],
};

// list of normal oxidation states for metals; these have quite a few exceptions, so they can only be considered guides
const OXSTATES:Record<string, number[]> =
{
	'Li': [1],
	'Na': [1],
	'K': [1],
	'Rb': [1],
	'Cs': [1],
	'Fr': [1],
	'Be': [2],
	'Mg': [2],
	'Ca': [2],
	'Sr': [2],
	'Ba': [2],
	'Ra': [2],
	'Al': [3],
	'Ga': [3],
	'In': [1, 3],
	'Tl': [1, 3],
	'Sn': [2, 4],
	'Pb': [2, 4],
	'Sb': [3, 5],
	'Bi': [3, 5],
	'Po': [2],
	'At': [1],
	'Sc': [3],
	'Y': [3],
	'Lu': [3],
	'Ti': [2, 4],
	'Zr': [4],
	'Hf': [4],
	'V': [5],
	'Nb': [5],
	'Ta': [5],
	'Cr': [2, 3, 6],
	'Mo': [0, 2, 3, 4, 6],
	'W': [2, 4, 6],
	'Mn': [1, 2, 3, 5, 7],
	'Tc': [7],
	'Re': [1, 3],
	'Fe': [0, 2, 3],
	'Ru': [0, 2, 3],
	'Os': [0, 2],
	'Co': [0, 1, 2, 3],
	'Rh': [1, 3],
	'Ir': [1, 3],
	'Ni': [0, 1, 2, 3],
	'Pd': [0, 2],
	'Pt': [0, 2],
	'Cu': [1, 2],
	'Ag': [1],
	'Au': [1, 2, 3],
	'Zn': [2],
	'Cd': [2],
	'Hg': [2],
};

export const enum AnalyseMoleculeType
{
	// problems
	BadValence = 'badvalence', // one of the 'boring' atoms has an invalid valence
	OddOxState = 'oxstate', // one of the 'metals' has an oxidation state that's suspicious
	WrongFormula = 'wrongformula', // calculated molecular formula doesn't match
	NonElement = 'nonelement', // not an atomic element

	// fixes
	FixCarbonyl = 'fixcarbonyl', // corrected a carbonyl that was missing unpaired
}

export interface AnalyseMoleculeResult
{
	type:AnalyseMoleculeType;
	atom?:number;
	bond?:number;
	value?:number;
	text?:string;
}

export class AnalyseMolecule
{
	public results:AnalyseMoleculeResult[] = [];
	public calcFormula:string;
	
	private atomMap:number[];
	private bondMap:number[];

	// ------------ public methods ------------

	// note that the molecule will be modified if necessary
	constructor(public mol:Molecule, private formula:string)
	{
	}

	public perform():AnalyseMoleculeResult[]
	{
		let mol = this.mol.clone();

		for (let n = 1; n <= mol.numAtoms; n++) mol.setAtomExtra(n, Vec.append(mol.atomExtra(n), 'xN:' + n));
		for (let n = 1; n <= mol.numBonds; n++) mol.setBondExtra(n, Vec.append(mol.bondExtra(n), 'xN:' + n));
		MolUtil.expandAbbrevs(mol, true);
		this.atomMap = Vec.numberArray(0, mol.numAtoms);
		this.bondMap = Vec.numberArray(0, mol.numBonds);
		for (let n = 1; n <= mol.numAtoms; n++) for (let ex of mol.atomExtra(n)) if (ex.startsWith('xN:')) this.atomMap[n - 1] = parseInt(ex.substring(3));
		for (let n = 1; n <= mol.numBonds; n++) for (let ex of mol.bondExtra(n)) if (ex.startsWith('xN:')) this.bondMap[n - 1] = parseInt(ex.substring(3));

		for (let n = 1; n <= mol.numAtoms; n++) this.identifyCorrections(mol, n);
		for (let n = 1; n <= mol.numAtoms; n++) this.identifyValenceProblems(mol, n);

		this.calcFormula = MolUtil.molecularFormula(mol);
		if (this.calcFormula != this.formula) this.results.push({'type': AnalyseMoleculeType.WrongFormula, 'text': this.calcFormula});

		return this.results;
	}

	// ------------ private methods ------------

	// find things that can be corrected
	private identifyCorrections(mol:Molecule, atom:number):void
	{
		let el = mol.atomElement(atom);
		if (el == 'C' && mol.atomUnpaired(atom) != 2 && mol.atomAdjCount(atom) == 2 && mol.atomHydrogens(atom) == 0)
		{
			let adj = mol.atomAdjBonds(atom);
			let o1 = mol.bondOrder(adj[0]), o2 = mol.bondOrder(adj[1]);
			if ((o1 == 0 && o2 == 2) || (o1 == 2 && o2 == 0))
			{
				mol.setAtomUnpaired(atom, 2);
				this.results.push({'type': AnalyseMoleculeType.FixCarbonyl, 'atom': this.atomMap[atom - 1]});
			}
		}
	}

	// find things that problematic (fatal or just noteworthy)
	private identifyValenceProblems(mol:Molecule, atom:number):void
	{
		if (mol.atomicNumber(atom) == 0)
		{
			this.results.push({'type': AnalyseMoleculeType.NonElement, 'atom': this.atomMap[atom - 1], 'text': mol.atomElement(atom)});
			return;;
		}

		let extra = mol.atomExtra(atom), extraMod = false;
		for (let i = extra.length - 1; i >= 0; i--) 
			if (extra[i].startsWith(BONDARTIFACT_EXTRA_RESPATH) || extra[i].startsWith(BONDARTIFACT_EXTRA_RESRING) || extra[i].startsWith(BONDARTIFACT_EXTRA_ARENE))
		{
			extra.splice(i, 1);
			extraMod = true;
		}
		if (extraMod) mol.setAtomExtra(atom, extra);

		if (mol.atomAdjCount(atom) == 0) return; // no point in making valence judgments for isolated atoms (?)

		let el = mol.atomElement(atom);
		let wantVal = VALENCES[el], wantOx = OXSTATES[el];
		if (!wantVal && !wantOx) return;

		// electronegativity ranking, for purposes of being "attached to an oxidant", to get the bonus higher oxidation levels; higher is more oxidising
		let eneg = (el:string):number =>
		{
			return ['Br', 'Cl', 'O', 'F'].indexOf(el);
		};

		let hyd = mol.atomHydrogens(atom), chg = mol.atomCharge(atom), unp = mol.atomUnpaired(atom);
		let valence = -chg + unp + hyd;
		let oxstate = chg + hyd;
		let oxiValence = 0; // record how much of the valence was due to O or F, which can trigger allowance of more states
		let tmBonds = Chemistry.ELEMENT_BLOCKS[mol.atomicNumber(atom)] >= 3 ? 0 : -1; // count # zero order bonds between transition metals
		let termDouble = 0; // count # of terminal double bonds (usually =O or =S) which usuallly count oxidation state differently
		for (let b of mol.atomAdjBonds(atom))
		{
			let order = mol.bondOrder(b), oatom = mol.bondOther(b, atom);
			valence += order;
			oxstate += order % 2;
			let otherEl = mol.atomElement(oatom);
			if (eneg(otherEl) > eneg(el)) oxiValence += order;
			if (tmBonds >= 0 && order == 0 && Chemistry.ELEMENT_BLOCKS[mol.atomicNumber(oatom)] >= 3) tmBonds++;
			if (order == 2 && mol.atomAdjCount(oatom) == 1 && mol.atomHydrogens(oatom) == 0 && mol.atomCharge(oatom) == 0) termDouble++;
		}

		// special deals on expanding valence options
		if (wantVal)
		{
			if (oxiValence > 0) wantVal = Vec.concat(wantVal, OXVALENCES[el]); // when oxidising neighbours
			if (chg > 0) wantVal = Vec.concat(wantVal, CATVALENCES[el]); // when it's a cation
		}
		if (wantVal && wantVal.indexOf(valence) < 0)
			this.results.push({'type': AnalyseMoleculeType.BadValence, 'atom': this.atomMap[atom - 1], 'value': valence});

		// oxidation state: ideally in the list; except that transition metal-metal bonds are conventionally written as zero-order, but
		// the mapping to common oxidation states often works better if we upgrade these to sigma bonds
		if (wantOx && wantOx.indexOf(oxstate) < 0 && wantOx.indexOf(oxstate + tmBonds) < 0 &&
					  wantOx.indexOf(oxstate + termDouble * 2) < 0 && wantOx.indexOf(oxstate + tmBonds + termDouble * 2) < 0)
			this.results.push({'type': AnalyseMoleculeType.OddOxState, 'atom': this.atomMap[atom - 1], 'value': oxstate});
	}
}

/* EOF */ }