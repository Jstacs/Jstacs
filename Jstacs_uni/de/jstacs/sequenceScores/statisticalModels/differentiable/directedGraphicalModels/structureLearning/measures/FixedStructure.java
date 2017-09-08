package de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures;

import de.jstacs.data.DataSet;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Class for a network structure of a
 * {@link de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM}
 * that has a given fixed dependency structure.
 * 
 * @author Jan Grau and Jens Keilwagen
 *
 */
public class FixedStructure extends Measure {
	
	//XXX problem: implements InstantiableFromParameterSet, but has no constructor

	private int[][] structure;

	/**
	 * Creates a new {@link FixedStructure} from its XML-representation.
	 * @param xml the XML-representation
	 * @throws NonParsableException the the XML could not be parsed
	 */
	public FixedStructure(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, getXMLTag() );
		structure = (int[][]) XMLParser.extractObjectForTags(xml, "structure");
	}

	/**
	 * The main constructor.
	 * 
	 * @param structure the dependency structure, <code>structure[i]</code> contains the parents of position <code>i</code>
	 */
	public FixedStructure(int[][] structure) {
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
		return getXMLTag();
	}

	@Override
	public int[][] getParents(DataSet fg, DataSet bg, double[] weightsFg, double[] weightsBg, int length)
			throws Exception {
		return ArrayHandler.clone(this.structure);
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, structure, "structure");
		XMLParser.addTags(xml, getXMLTag());
		return xml;
	}
}
