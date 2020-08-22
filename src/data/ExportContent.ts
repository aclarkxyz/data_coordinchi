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
///<reference path='../../../WebMolKit/src/data/MDLWriter.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Write the analysis molecules to an SDfile for use by other tools.
*/

export class ExportContent
{

	// ------------ public methods ------------

	constructor(private ds:DataSheet, private sdFN:string, private stereochemistry:boolean)
	{
	}

	public perform():void
	{
		console.log('Beginning export...');

		let colSame:number[] = [], colDiff:number[] = [];
		for (let n = 0; n < this.ds.numCols; n++) if (this.ds.colType(n) == DataSheetColumn.Molecule)
		{
			let cname = this.ds.colName(n);
			if (!cname.startsWith('Diff')) colSame.push(n); else colDiff.push(n);
		}

		let outDS = new DataSheet();
		let outMol = outDS.appendColumn('Molecule', DataSheetColumn.Molecule, '');
		let outGroup = outDS.appendColumn('Group', DataSheetColumn.Integer, '');

		let appendMolecule = (mol:Molecule, grp:number):void =>
		{
			let outRow = outDS.appendRow();
			outDS.setMolecule(outRow, outMol, mol);
			outDS.setInteger(outRow, outGroup, grp);

			// do permuted versions as well
			mol = mol.clone();
			for (let count = 3, n = 0; count > 0; count--, n++)
			{
				let a1 = (n % mol.numAtoms) + 1, a2 = ((n + 3) % mol.numAtoms) + 1;
				mol.swapAtoms(a1, a2);

				let outRow = outDS.appendRow();
				outDS.setMolecule(outRow, outMol, mol);
				outDS.setInteger(outRow, outGroup, grp);
			}
		};

		let grp = 0;
		for (let r = 0; r < this.ds.numRows; r++)
		{
			grp++;
			for (let c of colSame) if (this.ds.notNull(r, c))
			{
				let mol = this.ds.getMolecule(r, c);
				this.checkRoundTrip(mol, 'Row#' + (r + 1));
				appendMolecule(mol, grp);
			}

			for (let c of colDiff) if (this.ds.notNull(r, c))
			{
				let mol = this.ds.getMolecule(r, c);
				this.checkRoundTrip(mol, 'Row#' + (r + 1));
				appendMolecule(mol, ++grp);
			}
		}

		this.verifyAllHashes(outDS);

		let wtr = new MDLSDFWriter(outDS);
		wtr.enhancedFields = false;
		const fs = require('fs');
		fs.writeFileSync(this.sdFN, wtr.write());

		console.log('Done export.');
	}

	// ------------ private methods ------------

	// makes sure the molecule can survive an MDL Molfile round trip without using custom fields, which is necessary for
	// interoperability with other software
	private checkRoundTrip(mol:Molecule, src:string):void
	{
		mol = mol.clone();
		BondArtifact.removeAll(mol); // there is no molfile equivalent to this
		for (let n = 1; n <= mol.numBonds; n++)
		{
			let extra = mol.bondExtra(n).filter((str) => str != 'xAromatic' && str != 'xPi' && str != 'xDelocalised' && str != 'xUnknown');
			mol.setBondExtra(n, extra);
		}

		// pre-rotate abbreviations to the end, because parent-atom order isn't preserved
		let subidx = Vec.identity1(mol.numAtoms);
		for (let i = 0, limit = mol.numAtoms; i < limit;)
		{
			if (MolUtil.hasAbbrev(mol, subidx[i]))
			{
				subidx = Vec.remove(Vec.append(subidx, subidx[i]), i);
				limit--;
			}
			else i++;
		}
		mol = MolUtil.subgraphIndex(mol, subidx);

		let wtr = new MDLMOLWriter(mol);
		wtr.enhancedFields = false;
		wtr.abbrevSgroups = true;
		let molfile = wtr.write();
		let outmol = new MDLMOLReader(molfile).parse();

		let compareMols = (mol1:Molecule, mol2:Molecule):string =>
		{
			let na = mol1.numAtoms, nb = mol1.numBonds;
			if (na != mol2.numAtoms || nb != mol2.numBonds) return '# atoms';

			for (let n = 1; n <= na; n++)
			{
				if (mol1.atomElement(n) != mol2.atomElement(n)) return `element ${n}`;
				if (mol1.atomCharge(n) != mol2.atomCharge(n) || mol1.atomUnpaired(n) != mol2.atomUnpaired(n)) return `atom prop-a ${n}`;
				if (mol1.atomHydrogens(n) != mol2.atomHydrogens(n) || mol1.atomIsotope(n) != mol2.atomIsotope(n) || mol1.atomMapNum(n) != mol2.atomMapNum(n)) return `atom prop-b ${n}`;
				let tx1 = mol1.atomExtra(n).filter((x) => !x.startsWith('a'));
				let tx2 = mol2.atomExtra(n).filter((x) => !x.startsWith('a'));
				let ty1 = mol1.atomTransient(n), ty2 = mol2.atomTransient(n);
				tx1.sort(); tx2.sort(); ty1.sort(); ty2.sort();
				if (!Vec.equals(tx1, tx2) || !Vec.equals(ty1, ty2)) return `atom extra ${n}`;
			}

			let match = Vec.identity1(nb);
			for (let i = 1; i <= nb; i++)
			{
				let j = 0;
				for (let n = 0; n < match.length; n++)
				{
					j = match[n];
					if (mol1.bondFrom(i) == mol2.bondFrom(j) && mol1.bondTo(i) == mol2.bondTo(j)) {match.splice(n, 1); break;}
				}
				if (j == 0) return `bond ${i} unmatched`;

				let bo1 = Math.min(mol1.bondOrder(i), 3), bo2 = Math.min(mol2.bondOrder(j), 3);
				if (bo1 != bo2 || mol1.bondType(i) != mol2.bondType(j)) return `bond type ${i}/${j}`;
				let tx1 = mol1.bondExtra(i), tx2 = mol2.bondExtra(j), ty1 = mol1.bondTransient(i), ty2 = mol2.bondTransient(j);
				tx1.sort(); tx2.sort(); ty1.sort(); ty2.sort();
				if (!Vec.equals(tx1, tx2) || !Vec.equals(ty1, ty2)) return `bond extra ${i}/${j}`;
			}

			return null;
		};

		let diff = compareMols(mol, outmol);
		if (diff)
		{
			console.log('INPUT:\n' + mol);
			console.log('MOLFILE:\n' + molfile);
			console.log('OUTPUT:\n' + outmol);
			MolUtil.expandAbbrevs(mol, true);
			//console.log('EXPANDED ORIGINAL:\n' + mol);			
			throw `Read/write molfile failed [${diff}], source: ${src}`;
		}

		MolUtil.expandAbbrevs(mol, true);
		MolUtil.expandAbbrevs(outmol, true);
		let hash1 = new DotHash(new DotPath(mol), this.stereochemistry).calculate();
		let hash2 = new DotHash(new DotPath(outmol), this.stereochemistry).calculate();
		if (hash1 != hash2)
		{
			console.log('INPUT:\n' + mol);
			console.log('MOLFILE:\n' + molfile);
			console.log('OUTPUT:\n' + outmol);
			console.log('INHASH:  ' + hash1);
			console.log('OUTHASH: ' + hash2);
			throw 'Read/write molfile failed, dothashes different: ' + src;
		}
	}

	// final run-down before emitting the data: same hash iff same hash
	private verifyAllHashes(ds:DataSheet):void
	{
		let colMol = ds.findColByName('Molecule');
		let colGroup = ds.findColByName('Group');

		let molecules:Molecule[] = [], hashes:string[] = [];
		let sz = ds.numRows;
		for (let n = 0; n < sz; n++)
		{
			let mol = ds.getMoleculeClone(n, colMol);
			MolUtil.expandAbbrevs(mol, true);
			let hash = new DotHash(new DotPath(mol), this.stereochemistry).calculate();
			molecules.push(mol);
			hashes.push(hash);
		}

		for (let i = 0; i < sz - 1; i++) for (let j = i + 1; j < sz; j++)
		{
			let fail:string = null;
			if (ds.getInteger(i, colGroup) == ds.getInteger(j, colGroup))
			{
				if (hashes[i] != hashes[j]) fail = 'Hashes not matched: ' + [i + 1, j + 1];
			}
			else
			{
				if (hashes[i] == hashes[j]) fail = 'Hashes improperly matched: ' + [i + 1, j + 1];
			}
			if (fail)
			{
				console.log('Molecule 1:\n' + molecules[i]);
				console.log('Molecule 2:\n' + molecules[j]);
				console.log('Hash 1: ' + hashes[i]);
				console.log('Hash 2: ' + hashes[j]);
				throw fail;
			}
		}
	}
}

/* EOF */ }