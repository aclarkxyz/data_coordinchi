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
	DotHash: for a molecule that's represented by a DotPath, creates a unique hash code that can be used to distinguish between other equivalent
	structures.

	Assumptions: the provided DotPath instance assumes that the molecule has abbreviations expanded out. Having non-elements in the graph is allowed,
	but interpretation is undefined. Explicit/implicit/actual hydrogens do not affect the final hash code.
*/

export class DotHash
{
	public hash:string = null;

	private hcount:number[] = []; // # of implicit/explicit hydrogens (excluding actual hydrogens which have their own graph node)

	private atompri:number[]; // priority per-atom, which is resolved iteratively; final result: same priority value = completely identical
	private atomeqv:number[] = null; // atom equivalence: this is initially the same as atom priority, until "bumping" starts to disambiguate
	private chgNumer:number[]; // numerator for charges, divvied up into paths
	private chgDenom:number[]; // denominator of above
	private g:Graph; // neighbour-list
	private bondType:number[] = []; // bond-list analog: holds the dot-bond order
	private nbrType:number[][] = []; // neighbour-list analog: holds the dot-bond order

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
console.log('-----\n'+this.dot.mol);

		this.prepareMolecule();
		this.seedPriority();
		this.refinePriorities();
		this.composeHash();
//throw 'fnord';
		return this.hash;
	}

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

		let keepMask = Vec.booleanArray(true, mol.numAtoms);
		for (let n = 1; n <= na; n++) if (!pathMask[n - 1] && this.unwantedHydrogen(mol, n))
		{
			this.hcount[mol.atomAdjList(n)[0] - 1]++;
			keepMask[n - 1] = false;
		}
		if (Vec.anyFalse(keepMask))
		{
			mol = mol.clone();
			for (let n = 1; n <= na; n++) mol.setAtomHExplicit(n, this.hcount[n - 1]);

			mol = MolUtil.subgraphMask(mol, keepMask);
			this.hcount = Vec.maskGet(this.hcount, keepMask);
			this.dot = new DotPath(mol);
		}

		if (this.withStereo)
		{
			let meta = MetaMolecule.createRubric(mol);

			this.rubricTetra = Vec.anyArray(null, na);
			this.rubricSquare = Vec.anyArray(null, na);
			this.rubricBipy = Vec.anyArray(null, na);
			this.rubricOcta = Vec.anyArray(null, na);
			for (let n = 0; n < na; n++)
			{
				if (meta.rubricTetra[n]) this.rubricTetra[n] = Vec.sub(meta.rubricTetra[n], 1);
				else if (meta.rubricSquare[n]) this.rubricSquare[n] = Vec.sub(meta.rubricSquare[n], 1);
				else if (meta.rubricBipy[n]) this.rubricBipy[n] = Vec.sub(meta.rubricBipy[n], 1);
				else if (meta.rubricOcta[n]) this.rubricOcta[n] = Vec.sub(meta.rubricOcta[n], 1);
			}

			this.rubricSides = Vec.anyArray(null, na);
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

	// similar to the boringHydrogen method in MolUtil, except that it allows the elimination of stereoactive hydrogens; these are re-checked later; also allowed to
	// snip hydrogens bonded to exotic elements
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
		if (mol.atomHydrogens(other) + mol.atomAdjCount(other) > 4) return false; // neighbour has high valence

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

		/*console.log('PRISEQ:');
		for (let ps of priseq) console.log('    '+JSON.stringify(ps));
		console.log('PRI:'+JSON.stringify(this.atompri));
		console.log('NBR/order:'+JSON.stringify(this.nbrType));*/

	}

	// iteratively modify the atom priority array until uniqueness is achieved (by anti-symmetry bumping if necessary)
	private refinePriorities():void
	{
		const sz = this.atompri.length;
		while (true)
		{
			//let bits:string[] = [];
			//for (let n = 0; n < this.atompri.length; n++) bits.push(this.dot.mol.atomElement(n + 1) + ':' + this.atompri[n]);
			//console.log('ITER:'+JSON.stringify(bits));

			let adjpri:number[][] = [];
			for (let n = 0; n < sz; n++)
			{
				let rubric:number[] = null, perm:number[][] = null;
				if (!this.withStereo /*|| !this.atomeqv*/) {}
				else if (this.rubricTetra[n]) [rubric, perm] = [this.rubricTetra[n], Stereochemistry.RUBRIC_EQUIV_TETRA];
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

//console.log('---PREBUMP:'+this.atompri);
			this.bumpPriority();
//console.log('---POSTBUMP:'+this.atompri);
		}

		if (!this.atomeqv) this.defineEquivalents();
	}

	// when this is first called, the values of atompri indicate the canonical uniqueness of each atom, before taking into account stereochemistry; at this point it is
	// necessary to look through the "rubric" values and eliminate symmetrical cases
	private defineEquivalents():void
	{
		this.atomeqv = Vec.duplicate(this.atompri);
		if (!this.withStereo) return;

		const g = this.g, mol = this.dot.mol, na = mol.numAtoms, nb = mol.numBonds;

		for (let n = 0; n < na; n++)
		{
			// tetrahedral chirality inactive if any two equivalent neighbours
			if (this.rubricTetra[n])
			{
				let eq = this.neighbourWalkEquivalents([n], this.rubricTetra[n]);
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
				let eq = this.neighbourWalkEquivalents([n], this.rubricSquare[n]);
				if (multident) this.incorporateMultidentate(eq, n, this.rubricSquare[n]);
				for (let i = 0; i < 4; i++)
				{
					let nsame = 0;
					for (let j = 0; j < 4; j++) if (j != i && eq[i] == eq[j]) nsame++;
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
				let eq = this.neighbourWalkEquivalents([n], this.rubricBipy[n]);
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
				let eq = this.neighbourWalkEquivalents([n], this.rubricOcta[n]);
				if (multident) this.incorporateMultidentate(eq, n, this.rubricOcta[n]);
				for (let i = 0; i < 6; i++)
				{
					let nsame = 0;
					for (let j = 0; j < 6; j++) if (j != i && eq[i] == eq[j]) nsame++;
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
			let i1 = mol.bondFrom(n + 1) - 1, i2 = mol.bondTo(n + 1) - 1;
			let eq = this.neighbourWalkEquivalents([i1, i2], this.rubricSides[n]);
			if (eq[0] == eq[1] || eq[2] == eq[3]) this.rubricSides[n] = null;
		}
	}

	// for a list of neighbours immediately adjacent to atom(s), return a value for each such that same number = same thing; this takes into account meta
	// stereochemistry along the walk path
	private neighbourWalkEquivalents(src:number[], nbr:number[]):number[]
	{
		let pri:number[] = [];
		for (let n of nbr) pri.push(n < 0 ? 0 : this.atomeqv[n]);

/* (necessary?)
		const sz = nbr.length;
		let nsame = 0;
		for (let i = 0; i < sz - 1; i++) for (let j = i + 1; j < sz; j++) if (pri[i] == pri[j]) nsame++;
		if (nsame == 0) return pri; // they're all already different, no need to do any further analysis

		// examine each group
		let umask = Vec.booleanArray(false, sz);
		for (let i = 0; i < sz; i++) if (!umask[i])
		{
			let grp:number[] = [i];
			for (let j = i + 1; j < sz; j++) if (pri[j] == pri[j])
			{
				grp.push(j);
				umask[j] = true;
			}
		}

		// !! if they're all already different, stop
		// !! any group of same value, and not terminal, do a long walk out to see if they are really the same...
*/

		return pri;
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
			// ???? is this necessary ??
			//for (let n = 0; n < nbrpri.length; n++) nbrpri[n].unshift(n * 1000 + 1); // nbrpri[n].push(n + 1);
		}

		// NOTE: this would be where to insert the atom-centred stereochemistry constraint, to enforce an ordering
		let ret = [atompri[idx]];
		for (let pri of nbrpri) for (let p of pri) ret.push(p);
		return ret;
	}

	// find the lowest degenerate priority, and adjust it
	private bumpPriority():void
	{
		// part 1: consider stereo-actives first, and look for ways to disambiguate neighbours based on local geometry
		const sz = this.atompri.length;
		let pri = this.atompri.slice(0);
		for (let n = 0; n < sz; n++)
		{
			/*if (this.rubricTetra[n]) pri[n] += sz;
			else*/ if (this.rubricSquare[n]) pri[n] += 2 * sz;
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


console.log('considering:'+adjpri+'/'+Vec.uniqueUnstable(adjpri));
			if (Vec.uniqueUnstable(adjpri).length == adj.length) continue; // everything is already different, so move on

console.log('ATOM'+(n+1)+' atomRUBRIC:'+Vec.add(adj, 1));
			let extadj:number[][] = null;
			if (this.rubricSquare[n])
			{
				/*for (let i = 0; i < 4; i++)
				{
					let p1 = adjpri[i], p2 = adjpri[(i + 1) % 4], p3 = adjpri[(i + 2) % 4], p4 = adjpri[(i + 3) % 4];
					extadj.push([p1, Math.min(p2, p4), p3, Math.max(p2, p4)]);
				}*/

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
				/*let p1 = adjpri[0], p2 = adjpri[1], p3 = adjpri[2], p4 = adjpri[3], p5 = adjpri[4];
				extadj =
				[
					[p1, p2, p3, p4, p5],
					[p2, p3, p1, p4, p5],
					[p3, p1, p2, p4, p5],
					[p4, p5, p1, p2, p3],
					[p5, p4, p3, p2, p1],
				];
				const ROTEQUAT = [0, 2, 1, 4, 3];
				const ROTAXIAL = [[0, 1, 3, 4, 2], [0, 1, 4, 2, 3]];
				for (let i = 0; i < 3; i++)
				{
					let lookext = Vec.idxGet(extadj[i], ROTEQUAT);
					if (Vec.compareTo(lookext, extadj[i]) < 0) extadj[i] = lookext;
				}
				for (let i = 3; i < 5; i++)
				{
					let baseext = extadj[i];
					for (let perm of ROTAXIAL)
					{
						let lookext = Vec.idxGet(baseext, perm);
						if (Vec.compareTo(lookext, extadj[i]) < 0) extadj[i] = lookext;
					}
				}*/

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
				/*
				let p1 = adjpri[0], p2 = adjpri[1], p3 = adjpri[2], p4 = adjpri[3], p5 = adjpri[4], p6 = adjpri[5];
				extadj =
				[
					[p1, p2, p3, p4, p5, p6],
					[p2, p3, p4, p1, p5, p6],
					[p3, p4, p1, p2, p5, p6],
					[p4, p1, p2, p3, p5, p6],
					[p5, p1, p6, p3, p2, p4],
					[p6, p1, p5, p3, p4, p2],
				];
				const ROTATIONS = [[0, 3, 2, 1, 5, 4], [0, 4, 2, 5, 3, 1], [0, 5, 2, 4, 1, 3]];
				for (let i = 0; i < 6; i++)
				{
					let baseext = extadj[i];
					for (let perm of ROTATIONS)
					{
						let lookext = Vec.idxGet(baseext, perm);
						if (Vec.compareTo(lookext, extadj[i]) < 0) extadj[i] = lookext;
					}
				}*/

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
console.log('    EXTADJ:'+JSON.stringify(extadj));	//fnord
console.log('    OLDPRI:'+adjpri);
			let newpri = this.assignPriority(extadj);
console.log('    NEWPRI:'+newpri);

			let idx = Vec.idxSort(Vec.add(adjpri, Vec.mul(newpri, sz)));
console.log('    SORTED:'+idx+' NEWIDX:'+Vec.idxGet(newpri, idx));
			for (let i = 0; i < idx.length - 1; i++)
			{
				let i1 = idx[i], i2 = idx[i + 1];
				if (adj[i1] < 0) continue;
				if (adjpri[i1] != adjpri[i2] /*|| newpri[i1] == newpri[i2]*/) continue;

				let adjI = adj[i1];
console.log('       i='+i+' idx:'+[i1, i2]+' adjidx:'+adjI+','+adj[i2]);
				/*let bumpN = idx[i + 1];
				for (let j = 0; j < sz; n++) if (this.atompri
				this.atompri[bumpN]*/
				//this.atompri[idx[i]
console.log('      atompri_THEN:'+this.atompri);
				for (let j = 0; j < sz; j++) if (j != adjI && this.atompri[j] >= this.atompri[adjI]) this.atompri[j]++;
console.log('       atompri_NOW:'+this.atompri);
console.log('     rubricPRI:'+Vec.idxGet(this.atompri, adj));

				// priority has been purturbed, so push it around for another overall cycle
				return;
			}

//throw 'fnord';
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
			//if (debug) Util.writeln("Group Size=" + grp + "\nW" + (idx[n] + 1) + "=" + bestWalk.length + ":" + Util.arrayStr(bestWalk));

			for (let i = 1; i < grp; i++)
			{
				let walk = this.walkGraph(idx[n + i]);
				//if (debug) Util.writeln("W" + (idx[n + i] + 1) + "=" + walk.length + ":" + Util.arrayStr(walk));
				if (walk.length < bestWalk.length) {bestN = n + i; bestWalk = walk;}
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

			//if (debug) Util.writeln("Best=" + (idx[bestN] + 1));

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

				let stereoPri = this.stereoPriority(i2, pri, mask);

				shell.push([bondType[b - 1], pri[i2], stereoPri]);
				maskX[i2] = true;

				// !! rewrite below to order [i1,i2] to prev->next
				// !! in cases where rubric[next], permute the rubric so that mask[i] gets priority for the lower rungs, and !mask[i] goes next
				// !! incorporate the rubric into the shell sequence
				// (NOTE: if this works, check to see if incorporateMultident is still necessary...)

				/*if (mask[i1] && !mask[i2])
				{
					shell.push([bondType[b - 1], pri[i2]]);
					maskX[i2] = true;
				}
				else if (!mask[i1] && mask[i2])
				{
					shell.push([bondType[b - 1], pri[i1]]);
					maskX[i1] = true;
				}*/
			}
			for (let i = 0; i < sz; i++) mask[i] = maskX[i];

			this.sortSequences(shell);
			//for (let seq of shell) {path.push(seq[0]); path.push(seq[1]);}
			for (let seq of shell) path.push(seq[0], seq[1], seq[2]);
			path.push(-1);
		}
//console.log('PATHS:'+path);

		return path;
	}

	// for a given node index at some point along a breadth first walk, provide a differentiator based on stereo rubric: the mask indicates which of the atoms have
	// already been visited, and these are used as a priority
	private stereoPriority(idx:number, pri:number[], mask:boolean[]):number
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
	}

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

	// given that the atompri is defined and unique, composes a hash code that can disambiguate any two molecules
	private composeHash():void
	{
		const mol = this.dot.mol, na = mol.numAtoms, nb = mol.numBonds;

		let bits:string[] = [];

		// do atoms first
		let order = Vec.idxSort(this.atompri);
		for (let i of order)
		{
			let el = mol.atomElement(i + 1), hc = this.hcount[i], num = this.chgNumer[i], den = this.chgDenom[i];

			/*let par:string = null;
			if (!this.withStereo) {}
			else if (this.rubricTetra[i]) par = parity(-1, this.rubricTetra[i], Stereochemistry.RUBRIC_EQUIV_TETRA);
			else if (this.rubricSquare[i]) par = parity(i, this.rubricSquare[i], Stereochemistry.RUBRIC_EQUIV_SQUARE);
			else if (this.rubricBipy[i]) par = parity(i, this.rubricBipy[i], Stereochemistry.RUBRIC_EQUIV_BIPY);
			else if (this.rubricOcta[i]) par = parity(i, this.rubricOcta[i], Stereochemistry.RUBRIC_EQUIV_OCTA);
			else if (this.rubricSides[i]) par = parity(-1, this.rubricSides[i], [[0,1,2,3], [1,0,3,2]]);
			let pstr = par == null ? '' : '!' + par;*/
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

		this.hash = bits.join('');
	}

	// if the node index should have rubric parity associated with it, returns a string that can be inserted into the hash code as a disambiguation
	private parityString(idx:number):string
	{
		if (!this.withStereo) return null;
		/*let par:string = null;
		if (!this.withStereo) {}
		else if (this.rubricTetra[i]) par = parity(-1, this.rubricTetra[i], Stereochemistry.RUBRIC_EQUIV_TETRA);
		else if (this.rubricSquare[i]) par = parity(i, this.rubricSquare[i], Stereochemistry.RUBRIC_EQUIV_SQUARE);
		else if (this.rubricBipy[i]) par = parity(i, this.rubricBipy[i], Stereochemistry.RUBRIC_EQUIV_BIPY);
		else if (this.rubricOcta[i]) par = parity(i, this.rubricOcta[i], Stereochemistry.RUBRIC_EQUIV_OCTA);
		else if (this.rubricSides[i]) par = parity(-1, this.rubricSides[i], [[0,1,2,3], [1,0,3,2]]);
		let pstr = par == null ? '' : '!' + par;*/

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
			let g = this.g.clone(), mol = this.dot.mol;
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

			console.log('FROM:'+(idx+1)+' RUBRIC<atom>:'+Vec.add(rubric, 1));
			console.log('ADJ:'+JSON.stringify(optionAdj));
			console.log('PRI:'+JSON.stringify(optionPri));
			console.log('EQV:'+JSON.stringify(optionEqv));

			let order:number[]; 
			/*if (this.rubricBipy[idx])
			{
				order = [3, 4, 0, 1, 2];
				let lowpri = Math.min(optionPri[0][3], optionPri[0][4]);
				for (let n = optionAdj.length - 1; n >= 0; n--) if (optionPri[n][3] != lowpri)
				{
					optionAdj.splice(n, 1);
					optionPri.splice(n, 1);
					optionEqv.splice(n, 1);
				}
			}
			else
			{
				order = Vec.identity0(nsz);
				let lowpri = Vec.min(optionPri[0]);
				for (let n = optionAdj.length - 1; n >= 0; n--) if (optionPri[n][0] != lowpri)
				{
					optionAdj.splice(n, 1);
					optionPri.splice(n, 1);
					optionEqv.splice(n, 1);
				}
			}*/

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

			console.log('ORDER:'+order+' NPERM='+optionAdj.length);
			console.log('subADJ:'+JSON.stringify(optionAdj));
			console.log('subPRI:'+JSON.stringify(optionPri));
			console.log('subEQV:'+JSON.stringify(optionEqv));

			let blkgrp = new Map<number, number>(); // connected component -> group index
			let ngroups = 0;
			console.log('BLK:' + blk);
			if (blk[optionAdj[0][order[0]]] > 0) blkgrp.set(blk[optionAdj[0][order[0]]], ++ngroups);
			console.log('blkGroup:'+JSON.stringify(Array.from(blkgrp.entries())));

			let bestScore:number[];
			for (let n = 1; n < nsz; n++)
			{
	console.log('--- n='+n);		
				let bestN = -1, bestNode = -1;
				for (let i = 0; i < optionAdj.length; i++)
				{
					let score:number[] = [];
					for (let j = 0; j <= n; j++) 
					{
						let k = order[j], adj = optionAdj[i][k];
						let gscore = adj < 0 || blk[adj] == 0 ? 0 : (blkgrp.get(blk[adj]) || nsz);
	//console.log('        j='+j+' k='+k+' blk='+blk[optionAdj[i][k]] + ' gscore='+gscore);
						score.push(optionEqv[i][k] + sz * gscore);
					}
	console.log('    i='+i+' score='+score);				
					if (bestN < 0 || Vec.compareTo(score, bestScore) < 0) [bestN, bestScore, bestNode] = [i, score, optionAdj[i][order[n]]];
				}
	console.log('     bestNode='+bestNode+' blkBest='+blk[bestNode]+' blkgrphas='+blkgrp.has(blk[bestNode])+' ngr='+ngroups);
				if (blk[bestNode] > 0 && !blkgrp.has(blk[bestNode])) blkgrp.set(blk[bestNode], ++ngroups);
	console.log('     zblkGroup:'+JSON.stringify(Array.from(blkgrp.entries())));	
				for (let i = optionAdj.length - 1; i >= 0; i--) if (optionAdj[i][order[n]] != bestNode)
				{
					optionAdj.splice(i, 1);
					optionPri.splice(i, 1);
					optionEqv.splice(i, 1);
				}
			}

			console.log('outADJ:'+JSON.stringify(optionAdj));
			console.log('outPRI:'+JSON.stringify(optionPri));
			console.log('outEQV:'+JSON.stringify(optionEqv));

			console.log('bestScore:'+JSON.stringify(bestScore));

			return bestScore;
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

console.log('================= posFirst='+posFirst+' adjFirst='+adjFirst+'======');

		let bestScore:number[] = null;
		for (let n = 0; n < adjFirst.length; n++)
		{
console.log('n='+n+' adjFirst='+adjFirst[n]);		
			let score = generateParity(posFirst, adjFirst[n]);
			if (!score) continue;
console.log('        score='+score);
			if (bestScore == null || Vec.compareTo(score, bestScore) < 0) bestScore = score;
		}
console.log('BEST:'+bestScore);
/*/*if (this.rubricBipy[idx])
			{
				order = [3, 4, 0, 1, 2];
				let lowpri = Math.min(optionPri[0][3], optionPri[0][4]);
				for (let n = optionAdj.length - 1; n >= 0; n--) if (optionPri[n][3] != lowpri)
				{
					optionAdj.splice(n, 1);
					optionPri.splice(n, 1);
					optionEqv.splice(n, 1);
				}
			}
			else
			{
				order = Vec.identity0(nsz);
				let lowpri = Vec.min(optionPri[0]);
				for (let n = optionAdj.length - 1; n >= 0; n--) if (optionPri[n][0] != lowpri)
				{
					optionAdj.splice(n, 1);
					optionPri.splice(n, 1);
					optionEqv.splice(n, 1);
				}
			}*/


		let parity = Vec.idxSort(bestScore);
		console.log('parity:'+parity);
		//return parity;

		//throw 'fnord';

		return parity.toString();

		/*
		let adjpri = rubric.map((i) => rubric[i] < 0 ? 0 : this.atompri[rubric[i]]);
		let adjeqv = rubric.map((i) => rubric[i] < 0 ? 0 : this.atomeqv[rubric[i]]);

		// decide which neighbour to start with: normally just the lowest priority goes first, unless it's a lower symmetry configuration
		let pmap = Vec.numberArray(-1, nsz);
		if (this.rubricBipy[idx])
		{
			// pick one of the two axial positions
			if (adjpri[3] < adjpri[4])
				pmap[3] = 
			first = adjpri[3] < adjpri[4] ? 3 : 4;
		}
		else
		{
			for (let n = 0; n < nsz; n++) if (first < 0 || adjpri[n] < adjpri[first]) first = n;
			perms = perms.filter((perm) => perm[0] == first);
		}
*/

		//let nbrOrder:number[] = [first];
		//perms = perms.filter((perm) => perm[0] == first);

/*		let parity = (idx:number, nbr:number[], perm:number[][]):string =>
		{
			let adjpri = nbr.map((idx) => idx < 0 ? 0 : this.atomeqv[idx]);
			//let adjpri = nbr.map((idx) => idx < 0 ? 0 : this.atompri[idx]);

//console.log('atom='+(idx+1)+' adj='+adjpri+' alt='+nbr.map((idx) => idx < 0 ? 0 : this.atompri[idx])); // fnord
			//if (idx >= 0 && mol.atomRingBlock(idx + 1) > 0) this.incorporateMultidentate(adjpri, idx, nbr);
//console.log('  mult:'+adjpri);

			let bestpri = adjpri;
			for (let n = 1; n < perm.length; n++)
			{
				let permpri = Vec.idxGet(adjpri, perm[n]);
				if (Vec.compareTo(permpri, bestpri) < 0) bestpri = permpri;
			}

			// NOTE: returning the full parity array, sorted; most of the stereo types can be reduced to one bit of information; bipy & octa are a bit
			// more interesting though
			return Vec.idxSort(bestpri).toString();
			//return bestpri.toString();
		};*/

	}
}

/* EOF */ }