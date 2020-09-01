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

///<reference path='DotCompose.ts'/>

namespace WebMolKit /* BOF */ {

/*
	DotHash: for a molecule that's represented by a DotPath, creates a unique hash code that can be used to distinguish between other equivalent
	structures.

	Assumptions: the provided DotPath instance assumes that the molecule has abbreviations expanded out. Having non-elements in the graph is allowed,
	but interpretation is undefined. Explicit/implicit/actual hydrogens do not affect the final hash code.
*/

export class DotHash
{
	public hash:string = null;

	private hcount:number[] = []; // # of implicit/explicit hydrogens (excluding actual hydrogens which have their own graph node)
	private chgNumer:number[]; // numerator for charges, divvied up into paths
	private chgDenom:number[]; // denominator of above

	private atompri:number[]; // priority per-atom, which is resolved iteratively; final result: same priority value = completely identical
	private atomeqv:number[] = null; // atom equivalence: this is initially the same as atom priority, until "bumping" starts to disambiguate
	private g:Graph; // neighbour-list
	private bondType:number[] = []; // bond-list analog: holds the dot-bond order
	private nbrType:number[][] = []; // neighbour-list analog: holds the dot-bond order

	// zero-based rubric indices
	private rubricTetra:number[][] = null;
	private rubricSquare:number[][] = null;
	private rubricBipy:number[][] = null;
	private rubricOcta:number[][] = null;
	private rubricSides:number[][] = null;

	// ------------ public methods ------------

	// soft init: the hard work comes later
	constructor(public dot:DotPath, private withStereo:boolean)
	{
	}

	// build the hash, and return the calculated value
	public calculate():string
	{
		this.prepareMolecule();
		this.seedPriority();
		this.refinePriorities();

		let dotcomp = new DotCompose(this.dot, this.hcount, this.chgNumer, this.chgDenom, this.bondType, this.atompri, this.atomeqv);
		dotcomp.rubricTetra = this.rubricTetra;
		dotcomp.rubricSquare = this.rubricSquare;
		dotcomp.rubricBipy = this.rubricBipy;
		dotcomp.rubricOcta = this.rubricOcta;
		dotcomp.rubricSides = this.rubricSides;
		this.hash = dotcomp.compose();
		return this.hash;
	}

	public getAtomPrio():number[] {return this.atompri;}
	public getAtomEquiv():number[] {return this.atomeqv;}

	// ------------ private methods ------------

	// perform any necessary adjustments prior to creating the path
	private prepareMolecule():void
	{
		let mol = this.dot.mol, na = mol.numAtoms;

		// count up implicit hydrogens
		this.hcount = [];
		for (let n = 1; n <= mol.numAtoms; n++) this.hcount.push(mol.atomHydrogens(n));

		// see if any hydrogens need to be deleted
		let pathMask = Vec.booleanArray(false, mol.numAtoms);
		for (let pblk of this.dot.paths) for (let a of pblk.atoms) pathMask[a - 1] = true;

		let meta:MetaMolecule = null;
		if (this.withStereo) meta = MetaMolecule.createRubric(mol);

		let atomMask = Vec.booleanArray(true, mol.numAtoms), bondMask = Vec.booleanArray(true, mol.numBonds);
		for (let n = 1; n <= na; n++) if (!pathMask[n - 1] && this.unwantedHydrogen(mol, n))
		{
			if (meta)
			{
				// if the neighbour is involved in stereochemistry particular to heavy elements, keep it
				let idx = mol.atomAdjList(n)[0] - 1;
				if (meta.rubricSquare[idx] || meta.rubricBipy[idx] || meta.rubricOcta[idx]) continue;
			}

			this.hcount[mol.atomAdjList(n)[0] - 1]++;
			atomMask[n - 1] = false;
			bondMask[mol.atomAdjBonds(n)[0] - 1] = false;
		}
		if (Vec.anyFalse(atomMask))
		{
			mol = mol.clone();
			for (let n = 1; n <= na; n++) mol.setAtomHExplicit(n, this.hcount[n - 1]);

			mol = MolUtil.subgraphMask(mol, atomMask);
			this.hcount = Vec.maskGet(this.hcount, atomMask);
			this.dot = new DotPath(mol);

			// update the stereo-rubric too
			if (meta)
			{
				if (meta.rubricTetra) meta.rubricTetra = Vec.maskGet(meta.rubricTetra, atomMask);
				if (meta.rubricSquare) meta.rubricSquare = Vec.maskGet(meta.rubricSquare, atomMask);
				if (meta.rubricBipy) meta.rubricBipy = Vec.maskGet(meta.rubricBipy, atomMask);
				if (meta.rubricOcta) meta.rubricOcta = Vec.maskGet(meta.rubricOcta, atomMask);
				if (meta.rubricSides) meta.rubricSides = Vec.maskGet(meta.rubricSides, bondMask);

				let atomMap = Vec.prepend(Vec.add(Vec.maskMap(atomMask), 1), 0);
				for (let n = 0; n < Vec.arrayLength(meta.rubricTetra); n++) if (meta.rubricTetra[n]) meta.rubricTetra[n] = Vec.idxGet(atomMap, meta.rubricTetra[n]);
				for (let n = 0; n < Vec.arrayLength(meta.rubricSquare); n++) if (meta.rubricSquare[n]) meta.rubricSquare[n] = Vec.idxGet(atomMap, meta.rubricSquare[n]);
				for (let n = 0; n < Vec.arrayLength(meta.rubricOcta); n++) if (meta.rubricOcta[n]) meta.rubricOcta[n] = Vec.idxGet(atomMap, meta.rubricOcta[n]);
				for (let n = 0; n < Vec.arrayLength(meta.rubricSides); n++) if (meta.rubricSides[n]) meta.rubricSides[n] = Vec.idxGet(atomMap, meta.rubricSides[n]);
			}
		}

		this.rubricTetra = Vec.anyArray(null, na);
		this.rubricSquare = Vec.anyArray(null, na);
		this.rubricBipy = Vec.anyArray(null, na);
		this.rubricOcta = Vec.anyArray(null, na);
		this.rubricSides = Vec.anyArray(null, na);

		if (this.withStereo)
		{
			for (let n = 0; n < na; n++)
			{
				if (meta.rubricTetra[n]) this.rubricTetra[n] = Vec.sub(meta.rubricTetra[n], 1);
				else if (meta.rubricSquare[n]) this.rubricSquare[n] = Vec.sub(meta.rubricSquare[n], 1);
				else if (meta.rubricBipy[n]) this.rubricBipy[n] = Vec.sub(meta.rubricBipy[n], 1);
				else if (meta.rubricOcta[n]) this.rubricOcta[n] = Vec.sub(meta.rubricOcta[n], 1);
			}

			for (let n = 0; n < meta.rubricSides.length; n++)
			{
				let sides = meta.rubricSides[n];
				if (!sides || mol.bondInRing(n + 1)) continue;
				let [bfr, bto] = mol.bondFromTo(n + 1);
				this.rubricSides[bfr - 1] = [sides[0] - 1, sides[1] - 1, sides[2] - 1, sides[3] - 1];
				this.rubricSides[bto - 1] = [sides[2] - 1, sides[3] - 1, sides[0] - 1, sides[1] - 1];
			}
		}
	}

	// similar to the boringHydrogen method in MolUtil, except that it allows the elimination of stereoactive hydrogens; these are re-checked later; also disallowed from
	// snipping hydrogens bonded to exotic elements
	private unwantedHydrogen(mol:Molecule, atom:number):boolean
	{
		if (mol.atomElement(atom) != 'H') return false;

		if (mol.atomCharge(atom) != 0 || mol.atomUnpaired(atom) != 0) return false;
		if (mol.atomIsotope(atom) != Molecule.ISOTOPE_NATURAL) return false;
		if (mol.atomAdjCount(atom) != 1) return false;
		let other = mol.atomAdjList(atom)[0];
		if (mol.atomElement(other) == 'H') return false;
		let bond = mol.atomAdjBonds(atom)[0];
		if (mol.bondOrder(bond) != 1) return false;

		let atno = mol.atomicNumber(other);
		if (Chemistry.ELEMENT_BLOCKS[atno] >= 3) return false; // neighbour is d-block or later
		// ... this high-valence rule is not workable
		//if (mol.atomHydrogens(other) + mol.atomAdjCount(other) > 4) return false; // neighbour has high valence
		return true;
	}

	// setup the initial priority values prior to iterative convergence
	private seedPriority():void
	{
		const mol = this.dot.mol, na = mol.numAtoms, nb = mol.numBonds;

		// initial priorities for atoms: string together integer priorties in an array, then collapse it
		let priseq:number[][] = [];
		let chgNumer:number[] = [], chgDenom = Vec.numberArray(0, na);
		for (let n = 1; n <= na; n++) chgNumer.push(mol.atomCharge(n));
		for (let pblk of this.dot.paths)
		{
			let numer = 0, denom = pblk.atoms.length;
			for (let a of pblk.atoms) numer += mol.atomCharge(a);
			if (numer != 0)
			{
				let gcd = 1 / this.greatestCommonDenominator(Math.abs(numer), denom);
				numer *= gcd;
				denom *= gcd;
			}
			else denom = 0;

			for (let a of pblk.atoms) {chgNumer[a - 1] = numer; chgDenom[a - 1] = denom;}
		}
		for (let n = 1; n <= na; n++) priseq.push([mol.atomicNumber(n), this.hcount[n - 1], chgNumer[n - 1], chgDenom[n - 1]]);
		this.chgNumer = chgNumer;
		this.chgDenom = chgDenom;
		this.atompri = this.assignPriority(priseq);

		// bond orders: convert to neighbour list for quick access
		this.g = Graph.fromMolecule(mol);
		const g = this.g;
		this.bondType = this.dot.getBondClasses();
		for (let i = 0; i < na; i++)
		{
			let b:number[] = [];
			for (let j = 0; j < g.numEdges(i); j++)
			{
				let bidx = mol.findBond(i + 1, g.getEdge(i, j) + 1);
				b.push(this.bondType[bidx - 1]);
			}
			this.nbrType.push(b);
		}
	}

	// iteratively modify the atom priority array until uniqueness is achieved (by anti-symmetry bumping if necessary)
	private refinePriorities():void
	{
		const sz = this.atompri.length;
		while (true)
		{
			let adjpri:number[][] = [];
			for (let n = 0; n < sz; n++)
			{
				let rubric:number[] = null, perm:number[][] = null;
				if (this.rubricTetra[n]) [rubric, perm] = [this.rubricTetra[n], Stereochemistry.RUBRIC_EQUIV_TETRA];
				else if (this.rubricSquare[n]) [rubric, perm] = [this.rubricSquare[n], Stereochemistry.RUBRIC_EQUIV_SQUARE];
				else if (this.rubricBipy[n]) [rubric, perm] = [this.rubricBipy[n], Stereochemistry.RUBRIC_EQUIV_BIPY];
				else if (this.rubricOcta[n]) [rubric, perm] = [this.rubricOcta[n], Stereochemistry.RUBRIC_EQUIV_OCTA];
				else if (this.rubricSides[n]) [rubric, perm] = [this.rubricSides[n], [[0, 1, 2, 3], [1, 0, 3, 2]]];

				adjpri.push(this.adjacentPriority(n, rubric, perm));
			}

			let oldMax = Vec.max(this.atompri);
			this.atompri = this.assignPriority(adjpri);
			let newMax = Vec.max(this.atompri);

			if (newMax == this.atompri.length) break; // mission accomplished
			if (newMax > oldMax) continue; // achieved an incremental improvement

			// have to bump one of the priorities then let it percolate
			if (!this.atomeqv) this.defineEquivalents();

			this.bumpPriority();
		}

		if (!this.atomeqv) this.defineEquivalents();
	}

	// when this is first called, the values of atompri indicate the canonical uniqueness of each atom, before taking into account stereochemistry; at this point it is
	// necessary to look through the "rubric" values and eliminate symmetrical cases
	private defineEquivalents():void
	{
		this.atomeqv = Vec.duplicate(this.atompri);

		let neighbourWalkEquivalents= (nbr:number[]):number[] =>
		{
			let pri:number[] = [];
			for (let n of nbr) pri.push(n < 0 ? 0 : this.atomeqv[n]);
			return pri;
		};

		const g = this.g, mol = this.dot.mol, na = mol.numAtoms, nb = mol.numBonds;

		for (let n = 0; n < na; n++)
		{
			// tetrahedral chirality inactive if any two equivalent neighbours
			if (this.rubricTetra[n])
			{
				let eq = neighbourWalkEquivalents(this.rubricTetra[n]);
				outer: for (let i = 0; i < 3; i++) for (let j = i + 1; j < 4; j++) if (eq[i] == eq[j])
				{
					this.rubricTetra[n] = null;
					break outer;
				}
			}

			let multident = mol.atomRingBlock(n + 1) > 0; // things get more interesting when at least one ligand is multidentate

			// square planar stereochemistry inactive if any 3 identical constituents
			if (this.rubricSquare[n])
			{
				let eq = neighbourWalkEquivalents(this.rubricSquare[n]);
				if (multident) this.incorporateMultidentate(eq, n, this.rubricSquare[n]);
				for (let i = 0; i < 4; i++)
				{
					let nsame = 0;
					for (let j = 0; j < 4; j++) if (eq[i] == eq[j]) nsame++;
					if (nsame >= 3)
					{
						this.rubricSquare[n] = null;
						break;
					}
				}
			}

			/* this one is a bit different: if any of the 5 ligands are different, then there are multiple ways to assemble the ligand, so keep it
			// trigonal bipyramidal stereochemistry inactive if either of two orthonal planes
			if (this.rubricBipy[n])
			{
				let eq = this.neighbourWalkEquivalents([n], this.rubricBipy[n]);
				if (multident) this.incorporateMultidentate(eq, n, this.rubricBipy[n]);
				if (eq[3] == eq[4] || eq[0] == eq[3] || eq[0] == eq[2] || eq[1] == eq[2]) this.rubricBipy[n] = null;
			}*/
			// trigonal bipyramidal stereochemistry inactive only if all ligands are identical
			if (this.rubricBipy[n])
			{
				let eq = neighbourWalkEquivalents(this.rubricBipy[n]);
				if (multident) this.incorporateMultidentate(eq, n, this.rubricBipy[n]);
				let allSame = true;
				for (let i = 0; i < 4; i++) if (eq[i] != eq[i + 1])
				{
					allSame = false;
					break;
				}
				if (allSame) this.rubricBipy[n] = null;
			}

			// octahedral stereochemistry inactive if any 5 identical constituents
			if (this.rubricOcta[n])
			{
				let eq = neighbourWalkEquivalents(this.rubricOcta[n]);
				if (multident) this.incorporateMultidentate(eq, n, this.rubricOcta[n]);
				for (let i = 0; i < 6; i++)
				{
					let nsame = 0;
					for (let j = 0; j < 6; j++) if (eq[i] == eq[j]) nsame++;
					if (nsame >= 5)
					{
						this.rubricOcta[n] = null;
						break;
					}
				}
			}
		}

		// bondsides can be eliminated if either side is identical
		for (let n = 0; n < nb; n++) if (this.rubricSides[n])
		{
			//let i1 = mol.bondFrom(n + 1) - 1, i2 = mol.bondTo(n + 1) - 1;
			let eq = neighbourWalkEquivalents(this.rubricSides[n]);
			if (eq[0] == eq[1] || eq[2] == eq[3]) this.rubricSides[n] = null;
		}
	}

	// given that we have equivalence values for a stereocentre, factor in the possibility that we may have multidentate ligands that have the same
	// priorities, but by virtue of being part of a different ligand, they need to be distinguished; groups of unique multidentate ligands will get a
	// bump to their equivalence value which is arbitrary, but still preserves the canonical invariance of the calling context
	private incorporateMultidentate(eq:number[], centre:number, nbr:number[]):void
	{
		let g = Graph.fromMolecule(this.dot.mol);
		g.isolateNode(centre);
		let cc = g.calculateComponents();

		let ccToBlk = new Map<number, number[]>();
		for (let i of nbr) if (i >= 0)
		{
			let blk = ccToBlk.get(cc[i]);
			if (blk) blk.push(i); else ccToBlk.set(cc[i], [i]);
		}
		let groups:number[][] = [];
		for (let blk of ccToBlk.values()) if (blk.length > 1) groups.push(blk.map((i) => nbr.indexOf(i)));
		if (groups.length < 2) return; // no point

		let tags = groups.map((grp) => Vec.sorted(Vec.idxGet(eq, grp)).join(','));
		let bump = this.dot.mol.numAtoms;
		for (let i of Vec.idxSort(tags)) if (tags[i])
		{
			let idx = [i];
			for (let j = 0; j < groups.length; j++) if (i != j && tags[i] == tags[j]) idx.push(j);
			if (idx.length == 1) continue;
			for (let n of idx)
			{
				for (let j of groups[n]) eq[j] += bump;
				bump *= 2;
				tags[n] = null;
			}
		}
	}

	// for a 0-based node index, return an array that identifies its current priority sequence
	private adjacentPriority(idx:number, rubric:number[], perm:number[][]):number[]
	{
		const g = this.g, sz = g.numNodes, atompri = this.atompri, nbrType = this.nbrType;
		let nbrpri:number[][] = [];

		if (!rubric)
		{
			for (let n = 0; n < g.numEdges(idx); n++) nbrpri.push([nbrType[idx][n], atompri[g.getEdge(idx, n)]]);
			this.sortSequences(nbrpri);
		}
		else
		{
			for (let i of rubric)
			{
				if (i >= 0)
				{
					let j = g.getEdges(idx).indexOf(i);
					nbrpri.push([nbrType[idx][j], atompri[i]]);
				}
				else nbrpri.push([0, 0]);
			}
			this.sortPermutation(nbrpri, perm);
		}

		let ret = [atompri[idx]];
		for (let pri of nbrpri) for (let p of pri) ret.push(p);
		return ret;
	}

	// find the lowest degenerate priority, and adjust it
	private bumpPriority():void
	{
		// part 1: consider stereo-actives and look for ways to disambiguate neighbours based on local geometry
		const sz = this.atompri.length;
		let pri = this.atompri.slice(0);

		for (let n = 0; n < sz; n++)
		{
			//if (this.rubricTetra[n]) pri[n] += sz;
			//else
			if (this.rubricSquare[n]) pri[n] += 2 * sz;
			else if (this.rubricBipy[n]) pri[n] += 3 * sz;
			else if (this.rubricOcta[n]) pri[n] += 4 * sz;
			else pri[n] += 5 * sz; // non-stereoactives last
		}

		for (let n of Vec.idxSort(pri))
		{
			let adj:number[], perms:number[][];
			if (this.rubricSquare[n]) [adj, perms] = [this.rubricSquare[n], Stereochemistry.RUBRIC_EQUIV_SQUARE];
			else if (this.rubricBipy[n]) [adj, perms] = [this.rubricBipy[n], Stereochemistry.RUBRIC_EQUIV_BIPY];
			else if (this.rubricOcta[n]) [adj, perms] = [this.rubricOcta[n], Stereochemistry.RUBRIC_EQUIV_OCTA];
			else continue; // it's not one of the higher or stereotypes
			let adjpri = adj.map((a) => a < 0 ? 0 : this.atompri[a]);
			let [firstadj, firstpri] = [adj, adjpri];
			for (let perm of perms)
			{
				let lookadj = Vec.idxGet(firstadj, perm), lookpri = Vec.idxGet(firstpri, perm);
				if (Vec.compareTo(lookpri, adjpri) < 0) [adj, adjpri] = [lookadj, lookpri];
			}

			if (Vec.uniqueUnstable(adjpri).length == adj.length) continue; // everything is already different, so move on

			let extadj:number[][] = null;
			if (this.rubricSquare[n])
			{
				let bestPerm:number[] = null;
				for (let perm of Stereochemistry.RUBRIC_EQUIV_SQUARE)
				{
					let [p1, p2, p3, p4] = Vec.idxGet(adjpri, perm);
					let lookadj =
					[
						[p1, p2, p3, p4],
						[p2, p3, p4, p1],
						[p3, p4, p1, p2],
						[p4, p1, p2, p3],
					];
					if (extadj == null || Vec.compareTo(Vec.flatten(lookadj), Vec.flatten(extadj)) < 0) [extadj, bestPerm] = [lookadj, perm];
				}
				adj = Vec.idxGet(adj, bestPerm);
				adjpri = Vec.idxGet(adjpri, bestPerm);
			}
			else if (this.rubricBipy[n])
			{
				let bestPerm:number[] = null;
				for (let perm of Stereochemistry.RUBRIC_EQUIV_BIPY)
				{
					let [p1, p2, p3, p4, p5] = Vec.idxGet(adjpri, perm);
					let lookadj =
					[
						[p1, p2, p3, p4, p5],
						[p2, p3, p1, p4, p5],
						[p3, p1, p2, p4, p5],
						[p4, p5, p1, p2, p3],
						[p5, p4, p3, p2, p1],
					];
					if (extadj == null || Vec.compareTo(Vec.flatten(lookadj), Vec.flatten(extadj)) < 0) [extadj, bestPerm] = [lookadj, perm];
				}
				adj = Vec.idxGet(adj, bestPerm);
				adjpri = Vec.idxGet(adjpri, bestPerm);
			}
			else if (this.rubricOcta[n])
			{
				let bestPerm:number[] = null;
				for (let perm of Stereochemistry.RUBRIC_EQUIV_OCTA)
				{
					let [p1, p2, p3, p4, p5, p6] = Vec.idxGet(adjpri, perm);
					let lookadj =
					[
						[p1, p2, p3, p4, p5, p6],
						[p2, p3, p4, p1, p5, p6],
						[p3, p4, p1, p2, p5, p6],
						[p4, p1, p2, p3, p5, p6],
						[p5, p1, p6, p3, p2, p4],
						[p6, p1, p5, p3, p4, p2],
					];
					if (extadj == null || Vec.compareTo(Vec.flatten(lookadj), Vec.flatten(extadj)) < 0) [extadj, bestPerm] = [lookadj, perm];
				}
				adj = Vec.idxGet(adj, bestPerm);
				adjpri = Vec.idxGet(adjpri, bestPerm);
			}
			let newpri = this.assignPriority(extadj);

			let idx = Vec.idxSort(Vec.add(adjpri, Vec.mul(newpri, sz)));
			for (let i = 0; i < idx.length - 1; i++)
			{
				let i1 = idx[i], i2 = idx[i + 1];
				if (adj[i1] < 0) continue;
				if (adjpri[i1] != adjpri[i2]) continue;

				let adjI = adj[i1];
				for (let j = 0; j < sz; j++) if (j != adjI && this.atompri[j] >= this.atompri[adjI]) this.atompri[j]++;

				// priority has been perturbed, so push it around for another overall cycle
				return;
			}
		}

		// part 2: look for the canonically first occurring degenerates and perturb them
		pri = this.atompri;
		let idx = Vec.idxSort(pri);
		for (let n = 0; n < sz - 1; n++) if (pri[idx[n]] == pri[idx[n + 1]])
		{
			let grp = 2;
			while (n + grp < sz && pri[idx[n]] == pri[idx[n + grp]]) grp++;
			let bestN = n;
			let bestWalk = this.walkGraph(idx[n]);

			for (let i = 1; i < grp; i++)
			{
				let walk = this.walkGraph(idx[n + i]);
				if (walk.length < bestWalk.length) 
				{
					bestN = n + i;
					bestWalk = walk;
				}
				else if (walk.length == bestWalk.length)
				{
					for (let j = 0; j < walk.length; j++)
					{
						if (walk[j] < bestWalk[j])
						{
							bestN = n + i;
							bestWalk = walk;
							break;
						}
						if (walk[j] > bestWalk[j]) break;
					}
				}
			}
			for (; n < sz; n++) if (n != bestN) pri[idx[n]]++;
			break;
		}
	}

	// returns a sequence of indices that represent the BFS walk out through the graph
	private walkGraph(idx:number):number[]
	{
		const pri = this.atompri, sz = pri.length, mol = this.dot.mol, bondType = this.bondType
		let path:number[] = [];

		let mask = Vec.booleanArray(false, sz), maskX = Vec.booleanArray(false, sz);
		mask[idx] = true;

		for (let n = 0; n < sz; n++)
		{
			for (let i = 0; i < sz; i++) maskX[i] = mask[i];

			let shell:number[][] = [];
			for (let b = mol.numBonds; b >= 1; b--)
			{
				let i1 = mol.bondFrom(b) - 1, i2 = mol.bondTo(b) - 1;
				if (mask[i1] && !mask[i2]) {}
				else if (!mask[i1] && mask[i2]) [i1, i2] = [i2, i1];
				else continue;

				shell.push([bondType[b - 1], pri[i2], ...this.stereoNeighbours(i2)]);
				maskX[i2] = true;
			}
			for (let i = 0; i < sz; i++) mask[i] = maskX[i];

			this.sortSequences(shell);
			for (let seq of shell) path.push(...seq, -1);
			path.push(-2);
		}
		return path;
	}

	// returns the canonically smallest list of neighbour priorities that are compatible with the rubric, plus a type indicator: this can collapse
	// walk sequences that differ only by stereochemistry
	private stereoNeighbours(idx:number):number[]
	{
		let type:number, rubric:number[] = null, perms:number[][] = null;
		if (this.rubricTetra[idx]) [type, rubric, perms] = [1, this.rubricTetra[idx], Stereochemistry.RUBRIC_EQUIV_TETRA];
		else if (this.rubricSquare[idx]) [type, rubric, perms] = [2, this.rubricSquare[idx], Stereochemistry.RUBRIC_EQUIV_SQUARE];
		else if (this.rubricBipy[idx]) [type, rubric, perms] = [3, this.rubricBipy[idx], Stereochemistry.RUBRIC_EQUIV_BIPY];
		else if (this.rubricOcta[idx]) [type, rubric, perms] = [4, this.rubricOcta[idx], Stereochemistry.RUBRIC_EQUIV_OCTA];
		else return [];

		let sz = this.atompri.length;
		let rubricPri = rubric.map((i) => i < 0 ? 0 : this.atompri[i]);
		let bestPri:number[] = null;
		for (let perm of perms)
		{
			let pri = Vec.idxGet(rubricPri, perm);
			if (bestPri == null || Vec.compareTo(pri, bestPri) < 0) bestPri = pri;
		}
		return [type, ...bestPri];
	}

	// for a given node index at some point along a breadth first walk, provide a differentiator based on stereo rubric: the mask indicates which of the atoms have
	// already been visited, and these are used as a priority
	/*private stereoPriority(idx:number, pri:number[], mask:boolean[]):number
	{
		let rubric:number[] = null, perms:number[][] = null;
		if (this.rubricTetra[idx]) [rubric, perms] = [this.rubricTetra[idx], Stereochemistry.RUBRIC_EQUIV_TETRA];
		else if (this.rubricSquare[idx]) [rubric, perms] = [this.rubricSquare[idx], Stereochemistry.RUBRIC_EQUIV_SQUARE];
		else if (this.rubricBipy[idx]) [rubric, perms] = [this.rubricBipy[idx], Stereochemistry.RUBRIC_EQUIV_BIPY];
		else if (this.rubricOcta[idx]) [rubric, perms] = [this.rubricOcta[idx], Stereochemistry.RUBRIC_EQUIV_OCTA];
		else return 0;

		let bump = this.g.numNodes * 2;
		let basePri = rubric.map((r) => r < 0 ? 0 : pri[r] - (mask[r] ? bump : 0));
		let bestOrder = rubric, bestPri = basePri;

		skip: for (let n = 1; n < perms.length; n++)
		{
			let lookOrder = Vec.idxGet(rubric, perms[n]), lookPri = Vec.idxGet(basePri, perms[n]);
			for (let i = 0; i < basePri.length; i++)
			{
				if (lookPri[i] < bestPri[i]) break;
				if (lookPri[i] > bestPri[i]) continue skip;
			}
			[bestOrder, bestPri] = [lookOrder, lookPri];
		}

		let stereoPri = 0;
		bump = bestOrder.length;
		for (let n = 0; n < bestOrder.length; n++)
		{
			stereoPri += bestOrder[n] * bump;
			bump *= 2;
		}
		return stereoPri;
	}*/

	// for an array whereby each item is an array of numbers (arbitrary length), sorts these into a unique set and translates into an overall index, starting at 1
	private assignPriority(priseq:number[][]):number[]
	{
		let sorted = priseq.slice(0);
		this.sortSequences(sorted);
		let hashed:string[] = [];
		for (let seq of sorted)
		{
			let key = seq.toString();
			if (hashed.length == 0 || key != hashed[hashed.length - 1]) hashed.push(key);
		}
		let ret:number[] = [];
		for (let seq of priseq) ret.push(hashed.indexOf(seq.toString()) + 1);
		return ret;
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

	// select the lowest scoring permutation (equivalent to a constrained sort)
	private sortPermutation(priseq:number[][], perm:number[][]):void
	{
		let bestseq = priseq; // first permutation is identity
		for (let n = 1; n < perm.length; n++)
		{
			let permseq = Vec.idxGet(priseq, perm[n]);
			let cmp = 0;
			outer: for (let i = 0; i < priseq.length; i++)
			{
				let seq1 = bestseq[i], seq2 = permseq[i];
				const sz = Math.max(seq1.length, seq2.length);
				for (let j = 0; j < sz; j++)
				{
					let v1 = j < seq1.length ? seq1[j] : Number.NEGATIVE_INFINITY;
					let v2 = j < seq2.length ? seq2[j] : Number.NEGATIVE_INFINITY;
					if (v1 < v2) {cmp = -1; break outer;}
					if (v1 > v2) {cmp = 1; break outer;}
				}
			}
			if (cmp < 0) bestseq = permseq;
		}

		for (let n = 0; n < priseq.length; n++) priseq[n] = bestseq[n];
	}

	// quick & naive algorithm for finding the largest denominator of two numbers
	private greatestCommonDenominator(u:number, v:number):number
	{
		if (u < 0 || v < 0) throw 'Negative numbers disallowed.';
		if (u == v) return u;
		if (u == 0) return v;
		if (v == 0) return u;
		if (~u & 1) return (v & 1) ? this.greatestCommonDenominator(u >> 1, v) : this.greatestCommonDenominator(u >> 1, v >> 1) << 1;
		if (~v & 1) return this.greatestCommonDenominator(u, v >> 1);
		if (u > v) return this.greatestCommonDenominator((u - v) >> 1, v);
		return this.greatestCommonDenominator((v - u) >> 1, u);
	}


}

/* EOF */ }