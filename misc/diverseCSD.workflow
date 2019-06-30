<?xml version="1.0" encoding="UTF-8"?>
<Workflow>
	<Meta>
		<Description>Diverse CSD: selects a group of new candidates from CSD</Description>
		<Parameters>
			<Parameter name="infile" type="infile" descr="Input File" extensions="ds,sdf,el,mol,gz"/>
            <Parameter name="curfile" type="infile" descr="Current Selection" extensions="ds,sdf,el,mol,gz"/>
			<Parameter name="outfile" type="outfile" descr="Output File" extensions="ds"/>
            <Parameter name="divsize" type="integer" descr="Diverse subset to keep">100</Parameter>
		</Parameters>
	</Meta>

	<Nodes>
		<!-- read in the three inputs -->

		<Node name="ReadCSD" op="com.mmi.work.op.io.ReadAny">
			<Parameters>
				<Parameter name="filename" source="infile"/>
                <!--<Parameter name="maxRows">100000</Parameter>-->
			</Parameters>
			<Outputs count="1"/>
		</Node>

		<Node name="ReadExcl" op="com.mmi.work.op.io.ReadDataSheet">
			<Parameters>
				<Parameter name="filename" source="curfile"/>
			</Parameters>
			<Outputs count="1"/>
		</Node>

		<!-- first remove the current molecules -->

		<Node name="ExcludeCurrent" op="com.mmi.work.op.blk.DictionaryLookup">
            <Input name="ReadCSD" port="1"/>
            <Input name="ReadExcl" port="1"/>
            <Parameters>
                <Parameter name="keyLook">CSDID</Parameter>
                <Parameter name="keyFind">CSDID</Parameter>
                <Parameter name="excludeHit">true</Parameter>
            </Parameters>
            <Outputs count="1"/>
		</Node>

        <!-- pick a diverse subset -->

        <Node name="Diverse" op="com.mmi.work.op.mol.DiverseSubset">
            <Input name="ExcludeCurrent" port="1"/>
            <Parameters>
                <Parameter name="maxCount">100</Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>
		
        <!-- capture the results -->

		<Node name="WriteData" op="com.mmi.work.op.io.WriteDataSheet">
			<Parameters>
				<Parameter name="filename" source="outfile"/>
			</Parameters>
			<Input name="Diverse" port="1"/>
		</Node>
	</Nodes>
</Workflow>

