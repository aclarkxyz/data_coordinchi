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

		let bits:string[] = [];

		// do atoms first
		let order = Vec.idxSort(this.atompri);
		for (let i of order)
		{
			let el = mol.atomElement(i + 1), hc = this.hcount[i], num = this.chgNumer[i], den = this.chgDenom[i];
			let par = this.parityString(i);
			let pstr = par == null ? '' : '!' + par;

			bits.push('[' + el + ',' + hc + ',' + num + '/' + den + pstr + ']');
		}
		bits.push(';');

		// bonds next: they need to be sorted
		let backMap = Vec.numberArray(0, na);
		for (let n = 0; n < na; n++) backMap[order[n]] = n;
		let bondseq:number[][] = [];
		for (let n = 1; n <= nb; n++)
		{
			let a1 = backMap[mol.bondFrom(n) - 1] + 1;
			let a2 = backMap[mol.bondTo(n) - 1] + 1;
			[a1, a2] = [Math.min(a1, a2), Math.max(a1, a2)];
			bondseq.push([a1, a2, this.bondType[n - 1]]);
		}
		this.sortSequences(bondseq);
		for (let bs of bondseq) bits.push('[' + bs[0] + ':' + bs[1] + '=' + bs[2] + ']');

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
						//score.push(optionEqv[i][k] + sz * gscore);
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
		if (this.rubricTetra[idx] || this.rubricSides[idx])
		{
			return Permutation.parityIdentity(parity).toString();
		}
		else if (this.rubricSquare[idx])
		{
			// square pyramidal configurations can have 4 distinct permutation states, based on the scoring & ordering system
			if (Vec.equals(parity, [0, 1, 2, 3])) return '0';
			else if (Vec.equals(parity, [0, 1, 3, 2])) return '1';
			else if (Vec.equals(parity, [0, 2, 1, 3])) return '2';
			else if (Vec.equals(parity, [0, 2, 3, 1])) return '3';
			//else throw 'Invalid square planar parity: ' + parity + '/score=' + bestScore;
		}
		else if (this.rubricBipy[idx])
		{
			// trigonal bipyramidal can be distinguished by noting the priority indexing of the two axial ligands, followed by binary parity
			// indicating the order of the remaining substituents
			let invpar = [0, 0, 0, 0, 0];
			for (let n = 0; n < 5; n++) invpar[parity[n]] = n;
			return [invpar[3], invpar[4], Permutation.parityOrder([invpar[0], invpar[1], invpar[2]])].toString();
		}
		else if (this.rubricOcta[idx])
		{
			// octahedral can be distinguished by noting the two [arbitrary selected] "axial" ligands, followed by the binary parity
			// indicating the order of the remaining substituents
			let invpar = [0, 0, 0, 0, 0, 0];
			for (let n = 0; n < 6; n++) invpar[parity[n]] = n;
			return [invpar[4], invpar[5], Permutation.parityOrder([invpar[0], invpar[1], invpar[2], invpar[3]])].toString();
		}

		// NOTE: converting the parity array into a string works for all cases, but it includes more information than is necessary to
		// disambiguate the parities
		//return parity.toString();

		return null;
	}

	// sort in place: the array-of-arrays is sorted by lowest values first; does not assume that all elements are the same length
	private sortSequences(priseq:number[][]):void
	{
		priseq.sort((seq1, seq2):number =>
		{
			const sz = Math.max(seq1.length, seq2.length);
			for (let n = 0; n < sz; n++)
			{
				let v1 = n < seq1.length ? seq1[n] : Number.NEGATIVE_INFINITY;
				let v2 = n < seq2.length ? seq2[n] : Number.NEGATIVE_INFINITY;
				if (v1 < v2) return -1; else if (v1 > v2) return 1;
			}
			return 0;
		});
	}

}

/* EOF */ }