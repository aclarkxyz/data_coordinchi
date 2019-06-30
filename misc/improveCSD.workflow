<?xml version="1.0" encoding="UTF-8"?>
<Workflow>
	<Meta>
		<Description>Improve CSD: does some post-processing on imported structures</Description>
		<Parameters>
			<Parameter name="infile" type="infile" descr="Input File" extensions="ds,sdf,el,mol,gz"/>
			<Parameter name="outfile" type="outfile" descr="Output File" extensions="ds"/>
		</Parameters>
	</Meta>

	<Nodes>
		<Node name="ReadSelection" op="com.mmi.work.op.io.ReadDataSheet">
			<Parameters>
				<Parameter name="filename" source="infile"/>
			</Parameters>
			<Outputs count="1"/>
		</Node>

        <Node name="Improvement" op="com.mmi.db.pubinorg.InorgImprove">
            <Input name="ReadSelection" port="1"/>
            <Parameters>
            </Parameters>
            <Outputs count="1"/>
        </Node>

		<Node name="WriteImproved" op="com.mmi.work.op.io.WriteDataSheet">
			<Parameters>
				<Parameter name="filename" source="outfile"/>
			</Parameters>
			<Input name="Improvement" port="1"/>
		</Node>
	</Nodes>
</Workflow>

