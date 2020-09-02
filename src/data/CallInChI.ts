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
	private process:any;
	private listener:(inchi:string) => void = null;

	// ------------ public methods ------------

	constructor(private execPath:string, private withStereo:boolean)
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

	public async calculate(molList:Molecule[]):Promise<string[]>
	{
		if (!this.isAvailable) throw 'InChI unavailable';
		if (this.process) throw 'InChI is already in use';

		this.startProcess();

		return new Promise<string[]>((resolve, reject) =>
		{
			let inchiList:string[] = [];
			this.listener = (inchi:string) =>
			{
				inchiList.push(inchi);
				if (inchiList.length == molList.length)
				{
					this.process = null;
					resolve(inchiList);
				}
			};
			for (let mol of molList)
			{
				let feed = mol.clone();
				MolUtil.expandAbbrevs(feed, false);
				MolUtil.createHydrogens(feed, true);
				for (let n = 1; n <= feed.numBonds; n++) if (feed.bondOrder(n) == 0) feed.setBondOrder(n, 1);
				let wtr = new MDLMOLWriter(feed);
				wtr.enhancedFields = false;
				let mdl = wtr.write();
				this.process.stdin.write(mdl + '\n$$$$\n');
			}
			this.process.stdin.end();
		});
	}

	// ----------------- private methods -----------------

	private startProcess():void
	{
		const childproc = require('child_process'), path = require('path');

		let cmd = this.execPath.replace(/ /g, '\\\ '); // escape spaces
		let args = ['-STDIO', /*'-AuxNone',*/ '-NoLabels', /*'-Key',*/ '-DoNotAddH', '-FixedH', '-RecMet'];
		if (!this.withStereo) args.push('-SNon');
		let opt = {'stdio': ['pipe', 'pipe', 'pipe']};

		let lines:string[] = [];
		let buffer = '';

		this.process = childproc.spawn(cmd, args, opt);

		this.process.stdout.on('data', (chunk:Buffer | string) =>
		{
			buffer += chunk.toString();
			while (true)
			{
				let idx = buffer.indexOf('\n');
				if (idx < 0) break;
				lines.push(buffer.substring(0, idx));
				buffer = buffer.substring(idx + 1);
			}

			while (lines.length >= 2)
			{
				this.listener(lines[0]);
				lines.splice(0, 2);
			}
		});
		this.process.stderr.on('data', (chunk:Buffer | string) =>
		{
			//for (let line of (chunk.toString()).split('\n')) console.log('ERR: ' + line);
		});
	}
}

/* EOF */ }