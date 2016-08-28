package de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures;

import de.jstacs.data.DataSet;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;

public class FixedStructure extends Measure {

	private int[][] structure;

	public FixedStructure(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	public FixedStructure(int[][] structure) throws CloneNotSupportedException {
		super();
		this.structure = new int[structure.length][];
		for(int i=0;i<structure.length;i++){
			this.structure[i] = new int[structure[i].length+1];
			System.arraycopy(structure[i], 0, this.structure[i], 0, structure[i].length);
			this.structure[i][this.structure[i].length-1] = i;
		}
	}

	@Override
	public String getXMLTag() {
		return "Fixed";
	}

	@Override
	public String getInstanceName() {
		return "Fixed";
	}

	@Override
	public int[][] getParents(DataSet fg, DataSet bg, double[] weightsFg, double[] weightsBg, int length)
			throws Exception {
		return ArrayHandler.clone(this.structure);
	}

}
