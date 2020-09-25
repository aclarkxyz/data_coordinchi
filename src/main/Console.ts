/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/MDLReader.ts'/>

///<reference path='../decl/node.d.ts'/>

///<reference path='../data/DotHash.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Command line functionality.
*/

export class Console
{

	// ------------ public methods ------------

	constructor(private argv:string[])
	{
	}

	public async run():Promise<void>
	{
		if (this.argv.length < 1)
		{
			console.log('Reads an SDfile and outputs a list of coordination InChI layer identifiers to console');
			console.log('Command line:');
			console.log('    {input.sdf}');
			return;
		}

		let [ifn] = this.argv;
		let fs = require('fs');
		
		let ds:DataSheet;
		try {ds = new MDLSDFReader(fs.readFileSync(ifn).toString()).parse();}
		catch (ex)
		{
			console.log('Unable to read [' + ifn + ']: ' + ex);
			return;
		}

		let colMol = ds.firstColOfType(DataSheetColumn.Molecule);
		for (let n = 0; n < ds.numRows; n++)
		{
			let mol = ds.getMoleculeClone(n, colMol);
			MolUtil.expandAbbrevs(mol, true);
			let hash = new DotHash(new DotPath(mol), true).calculate();
			console.log(hash);
		}
	}

	// ------------ private methods ------------


}

// if NodeJS, needs to respond to exportion
if (typeof module == 'object') module.exports['Console'] = Console;

/* EOF */ }