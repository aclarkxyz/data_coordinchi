/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/sketcher/Sketcher.ts'/>
///<reference path='../../../WebMolKit/src/ui/ClipboardProxy.ts'/>
///<reference path='../../../WebMolKit/src/ui/Widget.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/MDLWriter.ts'/>
///<reference path='../../../WebMolKit/src/data/ForeignMolecule.ts'/>
///<reference path='../../../WebMolKit/src/gfx/Rendering.ts'/>
///<reference path='../../../WebMolKit/src/gfx/ArrangeMolecule.ts'/>
///<reference path='../../../WebMolKit/src/gfx/DrawMolecule.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

///<reference path='../data/AnalyseMolecule.ts'/>
///<reference path='../data/CallInChI.ts'/>
///<reference path='MoleculeCard.ts'/>
///<reference path='ZoomDotMol.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Equivalence results: reads a datasheet that has rows of structures that are supposed to be the same as each other, but not to any
	of the other rows.
*/

interface EquivalenceRow
{
	// molecule(s) expected to be the same
	molList:Molecule[];
	molExpanded:Molecule[];
	inchiList:string[];
	dhashList:string[];

	// molecule(s) expected to be different (in same row for clarity)
	diffList:Molecule[];
	diffExpanded:Molecule[];
	diffInChI:string[];
	diffDHash:string[];

	rubric:string;
}

export interface EquivalenceResultsOptions
{
	stereochemistry:boolean; // if true, stereoisomers are disambiguated
	failOnly:boolean; // only render rows that have a failure case of some kind
	inchiFail:boolean; // record failure when standard InChI fails to achieve desired effect (which happens a lot)
	startAt:number; // first row (1-based)
	endAt:number; // last row (1-based)
}

export class EquivalenceResults
{
	private cancelled = false;

	private divSummary:JQuery;
	private divContent:JQuery;
	private tableContent:JQuery;
	private tableRows:JQuery[] = [];

	private numPassed = 0;
	private numFailed = 0;
	private rowsFailed:number[] = [];
	private spanPassed:JQuery;
	private spanFailed:JQuery;
	private spanStatus:JQuery;
	private paraFailures:JQuery;

	private colMol:number[] = [];
	private colDiff:number[] = [];
	private colRubric = -1
	private rows:EquivalenceRow[] = [];

	protected policy:RenderPolicy;

	// ------------ public methods ------------

	constructor(private ds:DataSheet, private callInChI:CallInChI, private opt:EquivalenceResultsOptions, private callbackDone:() => void)
	{
	}

	public render(parentSummary:JQuery, parentContent:JQuery):void
	{
		this.divSummary = $('<div/>').appendTo(parentSummary);
		this.divContent = $('<div/>').appendTo(parentContent);

		this.policy = RenderPolicy.defaultColourOnWhite();
		this.policy.data.pointScale = 15;

		(async () => await this.startAnalysis())();
	}

	public cancel():void
	{
		this.cancelled = true;
	}

	// ------------ private methods ------------

	// get things setup, and start the chain of processing events
	private async startAnalysis():Promise<void>
	{
		// setup the summary panel
		let divStats = $('<div/>').appendTo(this.divSummary);
		divStats.append('Passed ');
		this.spanPassed = $('<span/>').appendTo(divStats).css({'border': '1px solid black', 'background-color': '#E0E0E0', 'padding': '0 0.25em 0 0.25em'});
		divStats.append(' Failed ');
		this.spanFailed = $('<span/>').appendTo(divStats).css({'border': '1px solid black', 'background-color': '#E0E0E0', 'padding': '0 0.25em 0 0.25em'});
		divStats.append(' ');
		this.spanStatus = $('<span/>').appendTo(divStats);
		this.paraFailures = $('<div/>').appendTo(this.divSummary);
		this.updateStats();

		// setup the table
		for (let n = 0; n < this.ds.numCols; n++) if (this.ds.colType(n) == DataSheetColumn.Molecule)
		{
			let cname = this.ds.colName(n);
			if (!cname.startsWith('Diff'))
				this.colMol.push(n);
			else
				this.colDiff.push(n);
		}
		this.colRubric = this.ds.findColByName('Rubric', DataSheetColumn.String);

		this.spanStatus.text('Preprocessing');

		let rosterMol:Molecule[] = [], rosterList:string[][] = [];

		for (let n = 0; n < this.ds.numRows; n++)
		{
			let eqr:EquivalenceRow =
			{
				'molList': [],
				'molExpanded': [],
				'inchiList': [],
				'dhashList': [],
				'diffList': [],
				'diffExpanded': [],
				'diffInChI': [],
				'diffDHash': [],
				'rubric': null,
			};
			for (let i = 0; i < this.colMol.length; i++) if (this.ds.notNull(n, this.colMol[i]))
			{
				let mol = this.ds.getMolecule(n, this.colMol[i]);
				let molExpanded = mol.clone();
				MolUtil.expandAbbrevs(molExpanded, true);

				// temporary bandaid: take the CSD-imported aromatic bond type and turn it into the "foreign" designation, which is used by dotpath
				for (let b = 1; b <= molExpanded.numBonds; b++) if (molExpanded.bondExtra(b).indexOf('xAromatic') >= 0)
				{
					let list = molExpanded.bondTransient(b);
					list.push(ForeignMoleculeExtra.BOND_AROMATIC);
					molExpanded.setBondTransient(b, list);
				}

				eqr.molList.push(mol);
				eqr.molExpanded.push(molExpanded);

				rosterMol.push(molExpanded);
				rosterList.push(eqr.inchiList);

				eqr.dhashList.push(null);
			}
			for (let i = 0; i < this.colDiff.length; i++) if (this.ds.notNull(n, this.colDiff[i]))
			{
				let mol = this.ds.getMolecule(n, this.colDiff[i]);
				let molExpanded = mol.clone();
				MolUtil.expandAbbrevs(molExpanded, true);

				eqr.diffList.push(mol);
				eqr.diffExpanded.push(molExpanded);

				rosterMol.push(molExpanded);
				rosterList.push(eqr.diffInChI);

				eqr.diffDHash.push(null);
			}
			if (this.colRubric >= 0) eqr.rubric = this.ds.getString(n, this.colRubric) || '';
			this.rows.push(eqr);
		}

		if (this.callInChI.isAvailable)
		{
			this.spanStatus.text('Calculating InChI');
			let rosterInChI = await this.callInChI.calculate(rosterMol);
			for (let n = 0; n < rosterInChI.length; n++) rosterList[n].push(rosterInChI[n]);
		}

		this.tableContent = $('<table/>').appendTo(this.divContent);

		let firstRow = this.opt.startAt - 1;
		if (!(firstRow >= 0)) firstRow = 0;

		this.processRow(firstRow);
	}

	private processRow(row:number):void
	{
		if (this.cancelled)
		{
			this.spanStatus.text('Cancelled');
			return;
		}

		if (row >= this.ds.numRows || (this.opt.endAt > 0 && row >= this.opt.endAt))
		{
			this.spanStatus.text('Finished');
			this.callbackDone();
			return;
		}

		let tr = $('<tr/>');
		this.tableRows[row] = tr;
		let th = $('<th/>').appendTo(tr).css({'text-align': 'left', 'vertical-align': 'top'});
		th.text('Row ' + (row + 1));

		let td = $('<td/>').appendTo(tr).css({'text-align': 'left', 'vertical-align': 'top', 'border': '1px solid #808080', 'background-color': '#F0F0F0'});
		let flex = $('<div/>').appendTo(td);
		flex.css({'display': 'flex', 'flex-wrap': 'wrap', 'justify-content': 'flex-start', 'align-items': 'flex-start'});

		let eqr = this.rows[row];

		if (eqr.rubric != null)
		{
			let strRubric = this.assembleRubric(eqr.molList[0]);
			if (strRubric != eqr.rubric)
			{
				console.log('** RUBRIC mismatch for row#' + (row + 1) + ' rubric:');
				console.log('    expected: ' + eqr.rubric);
				console.log('         got: ' + strRubric);
			}
		}

		let nmol = eqr.molList.length, ndiff = eqr.diffList.length;
		for (let n = 0; n < nmol; n++)
		{
			let mol = eqr.molList[n], molExpanded = eqr.molExpanded[n], inchi = eqr.inchiList[n];
			let dhash = new DotHash(new DotPath(molExpanded), this.opt.stereochemistry).calculate();
			eqr.dhashList[n] = dhash;

			let [card, spanMol] = this.generateCard(mol, molExpanded, inchi, dhash, 300);
			card.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #808080'});
			flex.append(card);
			spanMol.mouseenter(() => spanMol.css({'background-color': '#C0C0C0', 'border-radius': '5px'}));
			spanMol.mouseleave(() => spanMol.css('background-color', 'transparent'));
			spanMol.click(() => new ZoomDotMol(mol).open());

			// drill down on the hash codes in a bit more detail (sanity check)
			this.investigateHashes('Row ' + (row + 1) + '/Molecule ' + (n + 1), molExpanded, dhash);
		}
		for (let n = 0; n < ndiff; n++)
		{
			let molExpanded = eqr.diffExpanded[n];
			eqr.diffDHash[n] = new DotHash(new DotPath(molExpanded), this.opt.stereochemistry).calculate();
		}

		// compare equivalents within same row
		let hasProblem = false;
		for (let i = 0; i < nmol - 1; i++) for (let j = i + 1; j < nmol; j++)
		{
			let inchi1 = eqr.inchiList[i], inchi2 = eqr.inchiList[j];
			let dhash1 = eqr.dhashList[i], dhash2 = eqr.dhashList[j];
			let badInChI = inchi1 && inchi2 && inchi1 != inchi2, badHash = dhash1 != dhash2
			if (badInChI || badHash)
			{
				let mol1 = eqr.molList[i], mol2 = eqr.molList[j];
				let molExpanded1 = eqr.molExpanded[i], molExpanded2 = eqr.molExpanded[j];
				let card1 = this.generateCard(mol1, molExpanded1, inchi1, dhash1, 200)[0], card2 = this.generateCard(mol2, molExpanded2, inchi2, dhash2, 200)[0];

				let dualCard = $('<div/>').appendTo(flex);
				dualCard.css({'display': 'inline-block', 'margin': '0.5em'});
				dualCard.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #800000'});

				let divHdr = $('<div/>').appendTo(dualCard).css({'text-align': 'center', 'color': '#FF0000'});
				if (badInChI && badHash) divHdr.text('InChI & dots both different');
				else if (badInChI) divHdr.text('InChI codes differ');
				else if (badHash) divHdr.text('dot-hashes differ');

				let divMols = $('<div/>').appendTo(dualCard).css({'text-align': 'center'});
				divMols.append(card1);
				divMols.append(card2);
			}
			if (badHash || (badInChI && this.opt.inchiFail)) hasProblem = true;
		}

		// compare the reference molecule with any "differents" within the same row
		for (let n = 0; n < ndiff; n++)
		{
			let inchi1 = eqr.inchiList[0], inchi2 = eqr.diffInChI[n];
			let dhash1 = eqr.dhashList[0], dhash2 = eqr.diffDHash[n];
			let badInChI = inchi1 && inchi2 && inchi1 == inchi2, badHash = dhash1 == dhash2;
			if (badInChI || badHash)
			{
				let mol1 = eqr.molList[0], mol2 = eqr.diffList[n];
				let molExpanded1 = eqr.molExpanded[0], molExpanded2 = eqr.diffExpanded[n];
				let card1 = this.generateCard(mol1, molExpanded1, inchi1, dhash1, 200)[0], card2 = this.generateCard(mol2, molExpanded2, inchi2, dhash2, 200)[0];

				let dualCard = $('<div/>').appendTo(flex);
				dualCard.css({'display': 'inline-block', 'margin': '0.5em'});
				dualCard.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #800000'});

				let divHdr = $('<div/>').appendTo(dualCard).css({'text-align': 'center', 'color': '#FF0000'});
				if (badInChI && badHash) divHdr.text('InChI & dots both same');
				else if (badInChI) divHdr.text('InChI codes same');
				else if (badHash) divHdr.text('dot-hashes same');

				let divMols = $('<div/>').appendTo(dualCard).css({'text-align': 'center'});
				divMols.append(card1);
				divMols.append(card2);
			}
			if (badHash || (badInChI && this.opt.inchiFail)) hasProblem = true;
		}

		// compare between different rows: they are expected to be different
		for (let n = 0; n < this.rows.length; n++) if (n != row)
			for (let i = 0; i < nmol; i++) for (let j = 0; j < this.rows[n].molList.length; j++)
		{
			let other = this.rows[n];

			let inchi1 = eqr.inchiList[i], inchi2 = other.inchiList[j];
			let dhash1 = eqr.dhashList[i], dhash2 = other.dhashList[j];
			let badInChI = inchi1 && inchi2 && inchi1 == inchi2, badHash = dhash1 == dhash2

			if (badInChI || badHash)
			{
				let mol1 = eqr.molList[i], mol2 = other.molList[j];
				let molExpanded1 = eqr.molExpanded[i], molExpanded2 = other.molExpanded[j];
				let card1 = this.generateCard(mol1, molExpanded1, inchi1, dhash1, 200)[0], card2 = this.generateCard(mol2, molExpanded2, inchi2, dhash2, 200)[0];

				let dualCard = $('<div/>').appendTo(flex);
				dualCard.css({'display': 'inline-block', 'margin': '0.5em'});
				dualCard.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #800080'});

				let divHdr = $('<div/>').appendTo(dualCard).css({'text-align': 'center', 'color': '#800080'});
				let rowstr = ' [' + (row + 1) + ',' + (n + 1) + ']';
				if (badInChI && badHash) divHdr.text('InChI & dots falsely equivalent' + rowstr);
				else if (badInChI) divHdr.text('InChI codes falsely equivalent' + rowstr);
				else if (badHash) divHdr.text('Dot-hashes falsely equivalent' + rowstr);

				let divMols = $('<div/>').appendTo(dualCard).css({'text-align': 'center'});
				divMols.append(card1);
				divMols.append(card2);
			}
			if (badHash || (badInChI && this.opt.inchiFail)) hasProblem = true;
		}

		// by default add the row, unless the user has asked to skip failures
		if (hasProblem || !this.opt.failOnly) this.tableContent.append(tr);

		if (hasProblem)
		{
			this.numFailed++;
			this.rowsFailed.push(row);
		}
		else this.numPassed++;
		this.updateStats();

		setTimeout(() => this.processRow(row + 1), 1);
	}

	// try some variations on dotpath-hash generation; failure results in a hard crash, since this is a sanity check for the underlying algorithm, not a test of how
	// well the high level chemistry is working out
	private investigateHashes(note:string, mol:Molecule, dhash:string)
	{
		if (mol.numAtoms <= 1) return;

		let origMol = mol;
		mol = mol.clone();
		mol.keepTransient = true; // have some that we want to keep

		for (let count = 20, n = 0; count > 0; count--, n++)
		{
			let a1 = (n % mol.numAtoms) + 1, a2 = ((n + 3) % mol.numAtoms) + 1;
			mol.swapAtoms(a1, a2);
			let phash = new DotHash(new DotPath(mol), this.opt.stereochemistry).calculate();
			if (dhash != phash)
			{
				console.log('CONTENT:' + note + '/Iteration=' + (n + 1) + '/Perm=' + a1 + ':' + a2);
				console.log('ORIGINAL:' + dhash);
				console.log('PERMUTED:' + phash);
				console.log('Original Molecule:\n' + origMol);
				console.log('Permuted Molecule:\n' + mol);
				throw 'Dot hashes differ';
			}
		}
		//console.log('!! SKIP PERM');
	}

	// shorthand for creating a card object
	private generateCard(mol:Molecule, molExpanded:Molecule, inchi:string, dhash:string, dimsz:number):[JQuery, JQuery]
	{
		let card = new MoleculeCard(mol, molExpanded, inchi, dhash, dimsz, this.policy);
		card.generate();
		return [card.dom, card.spanMol];
	}

	// make sure panel is current
	private updateStats():void
	{
		this.spanPassed.text(this.numPassed.toString());
		this.spanFailed.text(this.numFailed.toString());
		this.spanStatus.text('Processing');

		this.paraFailures.empty();
		if (this.rowsFailed.length > 0)
		{
			this.paraFailures.append('Failure Rows:');
			for (let row of this.rowsFailed)
			{
				this.paraFailures.append(' ');
				let span = $('<span class="hover_action"/>').appendTo(this.paraFailures);
				span.text((row + 1).toString());
				span.click(() => $('html, body').animate({'scrollTop': this.tableRows[row].offset().top}, 500));
			}
		}
	}

	// produces a string that describes the stereochemistry rubric
	private assembleRubric(mol:Molecule):string
	{
		let bits:string[] = [];
		let meta = MetaMolecule.createRubric(mol);

		for (let n = 1; n <= mol.numAtoms; n++)
		{
			if (meta.rubricTetra[n - 1]) bits.push(`[T${n}:${meta.rubricTetra[n - 1]}]`);
			else if (meta.rubricSquare[n - 1]) bits.push(`[P${n}:${meta.rubricSquare[n - 1]}]`);
			else if (meta.rubricBipy[n - 1]) bits.push(`[B${n}:${meta.rubricBipy[n - 1]}]`);
			else if (meta.rubricOcta[n - 1]) bits.push(`[O${n}:${meta.rubricOcta[n - 1]}]`);
		}
		for (let n = 1; n <= mol.numBonds; n++) if (meta.rubricSides[n - 1] && !mol.bondInRing(n))
		{
			bits.push(`[S${n}:${meta.rubricSides[n - 1]}]`);
		}

		return bits.join(';');
	}
}

/* EOF */ }