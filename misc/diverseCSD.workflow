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

		<Node name="ReadCurrent" op="com.mmi.work.op.io.ReadDataSheet">
			<Parameters>
				<Parameter name="filename" source="curfile"/>
			</Parameters>
			<Outputs count="1"/>
		</Node>

        <!-- make multiple streams for the current content -->

        <Node name="SplitCurrent" op="com.mmi.work.op.blk.Multiplex">
            <Input name="ReadCurrent" port="1"/>
            <Outputs count="2"/>
        </Node>
        <Node name="Current1" op="com.mmi.work.op.blk.Buffer">
            <Input name="SplitCurrent" port="1"/>
            <Outputs count="1"/>
        </Node>
        <Node name="Current2" op="com.mmi.work.op.blk.Buffer">
            <Input name="SplitCurrent" port="2"/>
            <Outputs count="1"/>
        </Node>

        <!-- preliminary filtering -->

        <Node name="CalcProps" op="com.mmi.work.op.calc.MoleculeProperties">
            <Input name="ReadCSD" port="1"/>
            <Parameters>
                <Parameter name="heavyAtoms">HeavyAtoms</Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>

        <Node name="LimitHeavy" op="com.mmi.work.op.reduce.FilterProperties">
            <Input name="CalcProps" port="1"/>
            <Parameters>
                <Parameter name="column"><Item>HeavyAtoms</Item></Parameter>
                <Parameter name="operator"><Item>&lt;=</Item></Parameter>
                <Parameter name="value"><Item>100</Item></Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>

        <Node name="RemoveCols" op="com.mmi.work.op.fmt.DeleteColumns">
            <Input name="LimitHeavy" port="1"/>
            <Parameters>
                <Parameter name="columns"><Item>HeavyAtoms</Item></Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>

		<!-- first remove the current molecules -->

		<Node name="ExcludeCurrent" op="com.mmi.work.op.blk.DictionaryLookup">
            <Input name="RemoveCols" port="1"/>
            <Input name="Current1" port="1"/>
            <Parameters>
                <Parameter name="keyLook">CSDID</Parameter>
                <Parameter name="keyFind">CSDID</Parameter>
                <Parameter name="excludeHit">true</Parameter>
            </Parameters>
            <Outputs count="1"/>
		</Node>

        <!-- sort by similarity to current molecules -->

        <Node name="SimilarCurrent" op="com.mmi.work.op.mol.SimilarStructures">
            <Input name="ExcludeCurrent" port="1"/>
            <Input name="Current2" port="1"/>
            <Parameters>
                <Parameter name="column">Similarity</Parameter>
                <Parameter name="compare">maximum</Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>

        <Node name="SortSimilar" op="com.mmi.work.op.blk.Sort">
            <Input name="SimilarCurrent" port="1"/>
            <Parameters>
                <Parameter name="columns"><Item>Similarity</Item></Parameter>
                <Parameter name="directions"><Item>1</Item></Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>

        <Node name="KeepSimilar" op="com.mmi.work.op.reduce.KeepRows">
            <Input name="SortSimilar" port="1"/>
            <Parameters>
                <Parameter name="numKeep">100000</Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>

        <Node name="RemoveSimilar" op="com.mmi.work.op.fmt.DeleteColumns">
            <Input name="KeepSimilar" port="1"/>
            <Parameters>
                <Parameter name="columns"><Item>Similarity</Item></Parameter>
            </Parameters>
            <Outputs count="1"/>
        </Node>

        <!-- pick a diverse subset -->

        <Node name="Diverse" op="com.mmi.work.op.mol.DiverseSubset">
            <Input name="RemoveSimilar" port="1"/>
            <Parameters>
                <Parameter name="maxCount" source="divsize"/>
            </Parameters>
            <Outputs count="1"/>
        </Node>
		
        <!-- have a go at depicting -->

        <!-- ... insert the "improvement" feature ... -->

        <!-- capture the results -->

		<Node name="WriteData" op="com.mmi.work.op.io.WriteDataSheet">
			<Parameters>
				<Parameter name="filename" source="outfile"/>
			</Parameters>
			<Input name="Diverse" port="1"/>
		</Node>
	</Nodes>
</Workflow>

