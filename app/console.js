/*
	Coordination Complexes for InChI

	(c) 2019-2020 InChI Trust
*/

let WebMolKit = require(__dirname + '/coordinchi.js');
new WebMolKit.Console(process.argv.slice(2)).run().then();