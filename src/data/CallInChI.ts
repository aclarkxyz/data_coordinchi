/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/MolUtil.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Functionality for calling out to the InChI generator, which must be done by spawning a command line execution process.
*/

export interface CallInChIResult
{
	inchi:string;
	inchiKey:string;
}

export class CallInChI
{
	public isAvailable = false; // check this to see if it seems possible

	// ------------ public methods ------------

	constructor(private execPath:string)
	{
		if (!execPath) return;

		// see if the execution path corresponds to an executable file, and if not, fail silently
		const fs = require('fs');
		try
		{
			fs.accessSync(execPath, fs.constants.X_OK);
			this.isAvailable = true;
		}
		catch (ex) {} // not available
	}

	// fire up the executable and interpret the result; throws an exception if something went wrong
	public calculate(mol:Molecule):CallInChIResult
	{
		if (!this.isAvailable) throw 'The InChI generator is unconfigured or unavailable.';

		mol = mol.clone();
		MolUtil.expandAbbrevs(mol, false);
		MolUtil.createHydrogens(mol, true);

		const remote = require('electron').remote;
		const proc = require('child_process'), path = require('path');

		let cmd = this.execPath.replace(/ /g, '\\\ '); // escape spaces
		let writer = new MDLMOLWriter(mol);
		let mdlmol = writer.write();
		let result = proc.spawnSync(cmd, ['-STDIO', '-AuxNone', '-NoLabels', '-Key', '-DoNotAddH', '-SNon', '-FixedH', '-RecMet'], {'input': mdlmol});
		let raw = result.stdout.toString(), bits = raw.split('\n')

		if (bits.length < 2 || !bits[0].startsWith('InChI=1'))
		{
			console.log('** InChI generation failure');
			console.log('Output:\n' + raw);
			console.log('Input Molecule:\n' + mol);
			console.log('Input Molfile:\n' + mdlmol);
			throw 'InChI generation failed';
		}
		return {'inchi': bits[0], 'inchiKey': bits[1]};
	}

	// ----------------- private methods -----------------

}

/* EOF */ }