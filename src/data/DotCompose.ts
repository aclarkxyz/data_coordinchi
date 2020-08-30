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
///<reference path='../../../WebMolKit/src/data/Stereochemistry.ts'/>
///<reference path='../../../WebMolKit/src/data/DotPath.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Given that the DotHash class has worked out all of the important priority and ordering for a molecule, this class puts
	it together into a hashed string with InChI-esque characteristics.
*/

const BONDTYPE_MAP:[DotPathBond, string][] = 
[
	[DotPathBond.O0, '*'], [DotPathBond.O01, '*-'], [DotPathBond.O1, '-'], [DotPathBond.O12, '-='],
	[DotPathBond.O2, '='], [DotPathBond.O23, '=#'], [DotPathBond.O3, '#'], [DotPathBond.O3X, '#+']
];
const BONDTYPE_TEXT = new Map<DotPathBond, string>(BONDTYPE_MAP);

export class DotCompose
{
	public rubricTetra:number[][] = null;
	public rubricSquare:number[][] = null;
	public rubricBipy:number[][] = null;
	public rubricOcta:number[][] = null;
	public rubricSides:number[][] = null;	

	// ------------ public methods ------------

	constructor(private dot:DotPath, private hcount:number[], private chgNumer:number[], private chgDenom:number[], private bondType:number[],
				private atompri:number[], private atomeqv:number[])
	{
	}

	// build the hash, and return the calculated value
	public compose():string
	{
		const mol = this.dot.mol, na = mol.numAtoms, nb = mol.numBonds;
		const {bondType} = this;

		let bits:string[] = [];

		// do atoms first
		let atomLabel:string[] = [];
		for (let n = 0; n <na; n++)
		{
			let label = mol.atomElement(n + 1);
			if (this.hcount[n] > 0) label += 'H';
			if (this.hcount[n] > 1) label += this.hcount[n];
			if (this.chgNumer[n] > 0)
			{
				label += '+' + this.chgNumer[n];
				if (this.chgDenom[n] != 1) label += '/' + this.chgDenom[n];
			}
			let par = this.parityString(n);
			if (par) label += '!' + par;

			atomLabel.push(label);
		}

		let order = Vec.idxSort(this.atompri);

		for (let n = 0; n < order.length;)
		{
			if (bits.length > 0) bits.push(',');

			let blk = 1;
			for (; n + blk < order.length; blk++) if (atomLabel[order[n]] != atomLabel[order[n + blk]]) break;

			if (blk > 1) bits.push(blk + '*');
			bits.push(atomLabel[order[n]]);

			n += blk;
		}

		bits.push(';');

		// bonds next: they need to be sorted
		let backMap = Vec.numberArray(0, na); // graph index (0-based) to the ordered output atoms (1-based)
		for (let n = 0; n < na; n++) backMap[order[n]] = n + 1;

		let walkPaths = this.getWalkPaths();
		for (let n = 0; n < walkPaths.length; n++)
		{
			if (n > 0) bits.push(',');
			let path = walkPaths[n];
			bits.push(backMap[path[0] - 1].toString());
			for (let i = 1; i < path.length; i += 2)
			{
				let bond = path[i], atom = path[i + 1];
				bits.push(BONDTYPE_TEXT.get(bondType[bond - 1]));
				bits.push(backMap[atom - 1].toString());
			}
		}

		return bits.join('');
	}

	// ------------ private methods ------------

	// if the node index should have rubric parity associated with it, returns a string that can be inserted into the hash code as a disambiguation
	private parityString(idx:number):string
	{
		let rubric:number[], perms:number[][];
		if (this.rubricTetra[idx]) [rubric, perms] = [this.rubricTetra[idx], Stereochemistry.RUBRIC_EQUIV_TETRA];
		else if (this.rubricSquare[idx]) [rubric, perms] = [this.rubricSquare[idx], Stereochemistry.RUBRIC_EQUIV_SQUARE];
		else if (this.rubricBipy[idx]) [rubric, perms] = [this.rubricBipy[idx], Stereochemistry.RUBRIC_EQUIV_BIPY];
		else if (this.rubricOcta[idx]) [rubric, perms] = [this.rubricOcta[idx], Stereochemistry.RUBRIC_EQUIV_OCTA];
		else if (this.rubricSides[idx]) [rubric, perms] = [this.rubricSides[idx], [[0, 1, 2, 3], [1, 0, 3, 2]]];
		else return null;

		let sz = this.atompri.length, nsz = rubric.length;
		const {mol} = this.dot;
		let blk = Vec.numberArray(0, sz); // 0=not a bidentate ligand; >0=bidentate ligand identifier
		if (mol.atomRingBlock(idx + 1) > 0)
		{
			let mol = this.dot.mol, g = Graph.fromMolecule(mol);
			g.isolateNode(idx);
			let cc = g.calculateComponents();
			for (let adj of mol.atomAdjList(idx + 1)) if (mol.atomRingBlock(adj) == mol.atomRingBlock(idx + 1)) blk[adj - 1] = cc[adj - 1];
		}
		const {atompri, atomeqv} = this;

		let generateParity = (posFirst:number, adjFirst:number):number[] =>
		{
			let optionAdj:number[][] = [], optionPri:number[][] = [], optionEqv:number[][] = [];
			for (let perm of perms)
			{
				optionAdj.push(Vec.idxGet(rubric, perm));
				optionPri.push(perm.map((i) => rubric[i] < 0 ? 0 : atompri[rubric[i]]));
				optionEqv.push(perm.map((i) => rubric[i] < 0 ? 0 : atomeqv[rubric[i]]));
			}

			let order:number[];
			if (this.rubricBipy[idx])
				order = [3, 4, 0, 1, 2];
			else
				order = Vec.identity0(nsz);

			for (let n = optionAdj.length - 1; n >= 0; n--) if (optionAdj[n][posFirst] != adjFirst)
			{
				optionAdj.splice(n, 1);
				optionPri.splice(n, 1);
				optionEqv.splice(n, 1);
			}

			if (optionAdj.length == 0) return null;

			let blkgrp = new Map<number, number>(); // connected component -> group index
			let ngroups = 0;
			if (blk[optionAdj[0][order[0]]] > 0) blkgrp.set(blk[optionAdj[0][order[0]]], ++ngroups);

			let bestScore:number[];
			for (let n = 1; n < nsz; n++)
			{
				let bestN = -1, bestNode = -1;
				for (let i = 0; i < optionAdj.length; i++)
				{
					let score:number[] = [];
					for (let j = 0; j <= n; j++)
					{
						let k = order[j], adj = optionAdj[i][k];
						let gscore = adj < 0 || blk[adj] == 0 ? 0 : (blkgrp.get(blk[adj]) || nsz);
						score.push(sz * optionEqv[i][j] + gscore);
					}
					for (let j = 0; j <= n; j++) score.push(optionPri[i][order[j]]); // in case of a tie, prio kicks in
					if (bestN < 0 || Vec.compareTo(score, bestScore) < 0) [bestN, bestScore, bestNode] = [i, score, optionAdj[i][order[n]]];
				}
				if (blk[bestNode] > 0 && !blkgrp.has(blk[bestNode])) blkgrp.set(blk[bestNode], ++ngroups);
				for (let i = optionAdj.length - 1; i >= 0; i--) if (optionAdj[i][order[n]] != bestNode)
				{
					optionAdj.splice(i, 1);
					optionPri.splice(i, 1);
					optionEqv.splice(i, 1);
				}
			}

			return bestScore.slice(0, nsz);
		}

		// evaluate all equivalent starting points

		let posFirst = 0, adjFirst:number[] = [];
		if (this.rubricBipy[idx])
		{
			let loweqv = Math.min(atomeqv[rubric[3]], atomeqv[rubric[4]]);
			posFirst = 3;
			if (atomeqv[rubric[3]] == loweqv) adjFirst.push(rubric[3]);
			if (atomeqv[rubric[4]] == loweqv) adjFirst.push(rubric[4]);
		}
		else if (this.rubricSides[idx])
		{
			let loweqv = Number.POSITIVE_INFINITY;
			for (let n = 0; n < 2; n++) if (rubric[n] >= 0) loweqv = Math.min(loweqv, atomeqv[rubric[n]]);
			for (let n = 0; n < 2; n++) if (rubric[n] >= 0 && atomeqv[rubric[n]] == loweqv) adjFirst.push(rubric[n]);
		}
		else
		{
			let loweqv = Number.POSITIVE_INFINITY;
			for (let n = 0; n < nsz; n++) if (rubric[n] >= 0) loweqv = Math.min(loweqv, atomeqv[rubric[n]]);
			for (let n = 0; n < nsz; n++) if (rubric[n] >= 0 && atomeqv[rubric[n]] == loweqv) adjFirst.push(rubric[n]);
		}

		let bestScore:number[] = null;
		for (let n = 0; n < adjFirst.length; n++)
		{
			let score = generateParity(posFirst, adjFirst[n]);
			if (!score) continue;
			if (bestScore == null || Vec.compareTo(score, bestScore) < 0) bestScore = score;
		}

		let parity = Vec.idxSort(bestScore);

		// for tetrahedral or bond-side stereochemistry it is sufficient to represent the parity with a single bit; for the other stereoforms,
		// there are more degrees of freedom that need to be represented
		if (this.rubricTetra[idx])
		{
			return 't' + Permutation.parityIdentity(parity).toString();
		}
		else if (this.rubricSquare[idx])
		{
			// square pyramidal configurations can have 4 distinct permutation states, based on the scoring & ordering system
			if (Vec.equals(parity, [0, 1, 2, 3])) return 'q0';
			else if (Vec.equals(parity, [0, 1, 3, 2])) return 'q1';
			else if (Vec.equals(parity, [0, 2, 1, 3])) return 'q2';
			else if (Vec.equals(parity, [0, 2, 3, 1])) return 'q3';
			//else throw 'Invalid square planar parity: ' + parity + '/score=' + bestScore;
		}
		else if (this.rubricBipy[idx])
		{
			// trigonal bipyramidal can be distinguished by noting the priority indexing of the two axial ligands, followed by binary parity
			// indicating the order of the remaining substituents
			let invpar = [0, 0, 0, 0, 0];
			for (let n = 0; n < 5; n++) invpar[parity[n]] = n;
			return 'b' + [invpar[3], invpar[4], Permutation.parityOrder([invpar[0], invpar[1], invpar[2]])].join('');
		}
		else if (this.rubricOcta[idx])
		{
			// octahedral can be distinguished by noting the two [arbitrary selected] "axial" ligands, followed by the binary parity
			// indicating the order of the remaining substituents
			let invpar = [0, 0, 0, 0, 0, 0];
			for (let n = 0; n < 6; n++) invpar[parity[n]] = n;
			return 'o' + [invpar[4], invpar[5], Permutation.parityOrder([invpar[0], invpar[1], invpar[2], invpar[3]])].join('');
		}
		else if (this.rubricSides[idx])
		{
			return 's' + Permutation.parityIdentity(parity).toString();
		}

		// NOTE: converting the parity array into a string works for all cases, but it includes more information than is necessary to
		// disambiguate the parities
		//return parity.toString();

		return null;
	}

	// traverses the structure to generate some number of paths that traverses all of the bonds, in a canonical order, that tries to minimise
	// the number of paths required; the form of each path is [a1, b12, a2, b23, a3, ...] whereby each starts & ends with an atom (size = odd)
	private getWalkPaths():number[][]
	{
		const {atompri, bondType} = this;
		const {mol} = this.dot;
		const na = mol.numAtoms, nb = mol.numBonds;

		let bmask = Vec.booleanArray(false, nb); // set to true as each bond is visited and wrapped into a path
		let nbrList:number[][] = []; // zero-based atom adjacency which gets disconnected each time a bond is traversed
		for (let n = 1; n <= na; n++) nbrList.push(Vec.add(mol.atomAdjList(n), -1));

		// selects the next best atom for lowest walk path (0-based; -1 means it's all finished)
		let pickSeed = ():number =>
		{
			let nodes = Vec.identity0(na);

			// consider only atoms with at least one remaining adjacency
			nodes = nodes.filter((i) => nbrList[i].length > 0);
			if (nodes.length == 0) return -1; // the job is done

			// if any of these atoms is already part of a path, exclude all those which are not
			let usedMask = nodes.map((i) => nbrList[i].length < mol.atomAdjCount(i + 1));
			if (Vec.anyTrue(usedMask)) nodes = Vec.maskGet(nodes, usedMask);

			// find the minimum remaining adjacency count and limit to that
			let minadj = Vec.min(nodes.map((i) => nbrList[i].length));
			nodes = nodes.filter((i) => nbrList[i].length == minadj);

			// (could consider counting up the connected component size for each option and picking highest?)

			// of these remaining atoms, select the one with the lowest priority
			let lowidx = 0, lowpri = Number.POSITIVE_INFINITY;
			for (let i of nodes) if (atompri[i] < lowpri) [lowidx, lowpri] = [i, atompri[i]];
			return lowidx;
		};

		// considering the as-yet-unassigned bonds, and minus the current reference node, figures out how big the remaining connected
		// components are, which can bias the walk choices
		let compSize = (withoutAtom:number):number[] =>
		{
			let g = Graph.fromNeighbours(nbrList).clone();
			g.isolateNode(withoutAtom);
			let counts = Vec.numberArray(0, na);
			for (let cc of g.calculateComponentGroups()) for (let i of cc) counts[i] = cc.length;
			return counts;
		};

		let walkPaths:number[][] = [];

		for (let seed = pickSeed(); seed >= 0; seed = pickSeed())
		{
			let path = [seed]; // array is 0-based indexes while working on it, incremented to 1-based when it gets stashed

			while (true)
			{
				let head = Vec.last(path);
				let nbrs = nbrList[head];
				if (nbrs.length == 0) break;

				let next = nbrs[0];
				if (nbrs.length > 1)
				{
					// limit the available options to those which form the largest connected component
					let cmpsz = compSize(head);
					let maxsz = Vec.max(nbrs.map((i) => cmpsz[i]));
					nbrs = nbrs.filter((i) => cmpsz[i] == maxsz);
					next = nbrs[0];
					for (let n = 1; n < nbrs.length; n++) if (atompri[nbrs[n]] < atompri[next]) next = nbrs[n];
				}

				let bidx = mol.findBond(head + 1, next + 1) - 1;
				path.push(bidx);
				path.push(next);

				bmask[bidx] = true;
				nbrList[head] = nbrList[head].filter((i) => i != next);
				nbrList[next] = nbrList[next].filter((i) => i != head);
			}

			walkPaths.push(Vec.add(path, 1));
		}

		if (Vec.anyFalse(bmask)) throw 'Epic fail: ' + bmask;

		return walkPaths;
	}
}

/* EOF */ }