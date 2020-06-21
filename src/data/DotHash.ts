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
		this.prepareMolecule();
		this.seedPriority();
		this.refinePriorities();
		this.composeHash();

		return this.hash;
	}

	// ------------ private methods ------------

	// perform any necessary adjustments prior to creating the path
	private prepareMolecule():void
	{
		let mol = this.dot.mol, na = mol.numAtoms;

		if (this.withStereo)
		{
			let meta = MetaMolecule.createStrictRubric(mol);

			let zeroBased = (atoms:number[][]):number[][] =>
			{
				let zero = atoms.slice(0);
				for (let n = 0; n < zero.length; n++) if (zero[n]) zero[n] = Vec.sub(zero[n], 1);
				return zero;
			};

			this.rubricTetra = zeroBased(meta.rubricTetra);
			this.rubricSquare = zeroBased(meta.rubricSquare);
			this.rubricBipy = zeroBased(meta.rubricBipy);
			this.rubricOcta = zeroBased(meta.rubricOcta);

			this.rubricSides = Vec.anyArray(null, na);
			for (let n = 0; n < meta.rubricSides.length; n++)
			{
				let sides = meta.rubricSides[n];
				if (!sides) continue;
				let [bfr, bto] = mol.bondFromTo(n + 1);
				this.rubricSides[bfr - 1] = [sides[0] - 1, sides[1] - 1, sides[2] - 1, sides[3] - 1];
				this.rubricSides[bto - 1] = [sides[2] - 1, sides[3] - 1, sides[0] - 1, sides[1] - 1];
			}
		}

		// count up implicit hydrogens
		this.hcount = [];
		for (let n = 1; n <= mol.numAtoms; n++) this.hcount.push(mol.atomHydrogens(n));

		// see if any hydrogens need to be deleted
		let pathMask = Vec.booleanArray(false, mol.numAtoms);
		for (let pblk of this.dot.paths) for (let a of pblk.atoms) pathMask[a - 1] = true;

		let keepMask = Vec.booleanArray(true, mol.numAtoms);
		for (let n = 1; n <= na; n++) if (!pathMask[n - 1] && this.boringHydrogen(mol, n))
		{
			this.hcount[mol.atomAdjList(n)[0] - 1]++;
			keepMask[n - 1] = false;
		}
		if (Vec.anyFalse(keepMask))
		{
			mol = mol.clone(); // !! MAKE IT a bit faster
			for (let n = 1; n <= na; n++) mol.setAtomHExplicit(n, this.hcount[n - 1]);

			mol = MolUtil.subgraphMask(mol, keepMask);
			this.hcount = Vec.maskGet(this.hcount, keepMask);
			this.dot = new DotPath(mol);
		}
	}

	// same as the public method in MolUtil, except that it allows the elimination of stereoactive hydrogens; these are re-checked later; also allowed to
	// snip hydrogens bonded to exotic elements
	private boringHydrogen(mol:Molecule, atom:number):boolean
	{
		if (mol.atomElement(atom) != 'H') return false;

		if (mol.atomCharge(atom) != 0 || mol.atomUnpaired(atom) != 0) return false;
		if (mol.atomIsotope(atom) != Molecule.ISOTOPE_NATURAL) return false;
		if (mol.atomAdjCount(atom) != 1) return false;
		let other = mol.atomAdjList(atom)[0];
		if (mol.atomElement(other) == 'H') return false;
		let bond = mol.atomAdjBonds(atom)[0];
		if (mol.bondOrder(bond) != 1) return false;

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
		while (true)
		{
			//let bits:string[] = [];
			//for (let n = 0; n < this.atompri.length; n++) bits.push(this.dot.mol.atomElement(n + 1) + ':' + this.atompri[n]);
			//console.log('ITER:'+JSON.stringify(bits));

			let adjpri:number[][] = [];
			for (let n = 0; n < this.atompri.length; n++)
			{
				let rubric:number[] = null, perm:number[][] = null;
				if (!this.withStereo) {}
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
			this.bumpPriority();
		}

		if (!this.atomeqv) this.defineEquivalents();
	}

	// when this is first called, the values of atompri indicate the canonical uniqueness of each atom, before taking into account stereochemistry; at this point it is
	// necessary to look through the "rubric" values and eliminate symmetrical cases
	private defineEquivalents():void
	{
		this.atomeqv = Vec.duplicate(this.atompri);
		if (!this.withStereo) return;

		// !! TODO: look for symmetry within the rubric, to eliminate some o fthem
	}

	// for a 0-based node index, return an array that identifies its current priority sequence; if the rubric/perm parameters are not null, these
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
//console.log('ADJACENT:atom='+(idx+1)+' rubric='+rubric); //fnord
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
		const pri = this.atompri, sz = pri.length;
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
				if (mask[i1] && !mask[i2])
				{
					shell.push([bondType[b - 1], pri[i2]]);
					maskX[i2] = true;
				}
				else if (!mask[i1] && mask[i2])
				{
					shell.push([bondType[b - 1], pri[i1]]);
					maskX[i1] = true;
				}
			}
			for (let i = 0; i < sz; i++) mask[i] = maskX[i];

			this.sortSequences(shell);
			for (let seq of shell) {path.push(seq[0]); path.push(seq[1]);}
			path.push(0);
		}

		return path;
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
//console.log('SORT:'+JSON.stringify(priseq));
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
//console.log(' CMP:'+JSON.stringify(permseq)+' cmp='+cmp);
			if (cmp < 0) bestseq = permseq;
		}

//fnord
//console.log('BEST:'+JSON.stringify(bestseq)+'\n');

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

		/*let parity = (nbr:number[], perm:number[][]):string =>
		{
			let adjpri = nbr.map((idx) => idx < 0 ? 0 : this.atomeqv[idx]);

			let bestpri = adjpri;
			for (let n = 1; n < perm.length; n++)
			{
				let permpri = Vec.idxGet(adjpri, perm[n]);
				for (let i = 0; i < permpri.length; i++) 
				{
					if (permpri[i] < bestpri[i]) {bestpri = permpri; break;}
					if (permpri[i] > bestpri[i]) break;
				}
			}

			let p = 0;
			for (let i = 0; i < nbr.length - 1; i++) for (let j = i + 1; j < nbr.length; j++) if (bestpri[i] > bestpri[j]) p++;
			return (p & 1).toString();
		};*/
		let parity = (nbr:number[], perm:number[][]):string =>
		{
			let adjpri = nbr.map((idx) => idx < 0 ? 0 : this.atomeqv[idx]);

			let bestpri = adjpri;
			for (let n = 1; n < perm.length; n++)
			{
				let permpri = Vec.idxGet(adjpri, perm[n]);
				for (let i = 0; i < permpri.length; i++) 
				{
					if (permpri[i] < bestpri[i]) {bestpri = permpri; break;}
					if (permpri[i] > bestpri[i]) break;
				}
			}

			// NOTE: returning the full parity array, sorted; most of the stereo types can be reduced to one bit of information; bipy & octa are a bit
			// more interesting though
			return Vec.idxSort(bestpri).toString();
		};

		let bits:string[] = [];

//console.log('PRIORITIES:' + this.atompri);
//console.log('EQUIVS:'+this.atomeqv);
		// do atoms first
		let order = Vec.idxSort(this.atompri);
//console.log('-----COMPOSING----')		
		for (let i of order)
		{
			let el = mol.atomElement(i + 1), hc = this.hcount[i], num = this.chgNumer[i], den = this.chgDenom[i];

			let par:string = null;
//if (this.rubricTetra[i]) console.log('ATOM:'+i+' pri='+this.atompri[i]+' eqv='+this.atomeqv[i]+' adj='+this.g.getEdges(i)+' tetra='+this.rubricTetra[i]);//fnord
			if (!this.withStereo) {}
			else if (this.rubricTetra[i]) par = parity(this.rubricTetra[i], Stereochemistry.RUBRIC_EQUIV_TETRA);
			else if (this.rubricSquare[i]) par = parity(this.rubricSquare[i], Stereochemistry.RUBRIC_EQUIV_SQUARE);
			else if (this.rubricBipy[i]) par = parity(this.rubricBipy[i], Stereochemistry.RUBRIC_EQUIV_BIPY);
			else if (this.rubricOcta[i]) par = parity(this.rubricOcta[i], Stereochemistry.RUBRIC_EQUIV_OCTA);
			else if (this.rubricSides[i]) par = parity(this.rubricSides[i], [[0,1,2,3], [1,0,3,2]]);			
			let pstr = par == null ? '' : '!' + par;
//console.log('atom='+(i+1)+'/'+mol.atomElement(i + 1)+' '+JSON.stringify([this.rubricTetra[i], this.rubricSquare[i], this.rubricBipy[i], this.rubricOcta[i], this.rubricSides[i]])+ ' parity='+par);

			bits.push('[' + el + ',' + hc + ',' + num + '/' + den + pstr + ']');
		}
//throw 'fnord'; // zog
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
}

/* EOF */ }