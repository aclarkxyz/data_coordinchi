/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/DataSheetStream.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>
///<reference path='../data/AnalyseMolecule.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Command line: analyze a datasheet, by going through and looking for problems with each of the structures.
*/

export class AnalyseDataSheets
{
	private filenames:string[] = [];

	// ------------ public methods ------------

	constructor(args:string[])
	{
		for (let n = 0; n < args.length; n++)
		{
			if (args[n].startsWith('-')) throw 'Unexpected argument: ' + args[n];
			else this.filenames.push(args[n]);
		}
		if (this.filenames.length == 0) throw 'Must provide at least one datasheet filename argument.';
	}

	public exec():void
	{
		console.log('Analyzing files...');

		const fs = require('fs');

		for (let n = 0; n < this.filenames.length; n++)
		{
			let fn = this.filenames[n];
			console.log('  ' + (n + 1) + '/' + this.filenames.length + ': ' + fn);
			let strXML = '';
			try {strXML = fs.readFileSync(fn);}
			catch (ex) {throw 'Unable to read file: ' + fn;}
			let ds = DataSheetStream.readXML(fn);
			if (ds == null) throw 'Unable to parse file ' + fn;

			this.analyse(ds);

			// !! write...
		}

		console.log('Done.');
	}

	// ------------ private methods ------------

	private analyse(ds:DataSheet):void
	{
		let colMol = ds.firstColOfType(DataSheet.COLTYPE_MOLECULE);
		let colError = ds.ensureColumn('Errors', DataSheet.COLTYPE_STRING, 'Fatal flaws with the incoming molecule');
		let colWarning = ds.ensureColumn('Warnings', DataSheet.COLTYPE_STRING, 'Questionable attributes of the incoming molecule');

		for (let n = 0; n < ds.numRows; n++)
		{
			let mol = ds.getMolecule(n, colMol);
			let anal = new AnalyseMolecule(mol);
			anal.perform();

			let error:string[] = [], warning:string[] = [];
			for (let result of anal.results)
			{
				if (result.type == AnalyseMoleculeType.BadValence)
				{
					error.push('Valence:atom=' + result.atom + '[' + mol.atomElement(result.atom) + ']:' + result.value);
				}
				else if (result.type == AnalyseMoleculeType.OddOxState)
				{
					warning.push('OxState:atom=' + result.atom + '[' + mol.atomElement(result.atom) + ']:' + result.value);
				}
			}

			ds.setMolecule(n, colMol, anal.mol);
			ds.setString(n, colError, error.join('\n'));
			ds.setString(n, colWarning, warning.join('\n'));
		}
	}
}

/* EOF */ }