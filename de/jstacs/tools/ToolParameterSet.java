package de.jstacs.tools;

import java.util.Collection;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.ResultSetResult;

/**
 * Class for the parameters of a {@link JstacsTool}.
 * In addition to a standard {@link ResultSetResult}, this class stores some supplementary information
 * like the name of the tool.
 *  
 * @author Jens Keilwagen, Jan Grau
 * 
 * @see JstacsTool
 * @see JstacsTool#getToolParameters()
 */
public class ToolParameterSet extends ParameterSet {
	
	protected String toolName;
	
	/**
	 * Constructs a {@link ToolParameterSet} given a tool name and some {@link Parameter}s.
	 * The {@link Parameter}s are not cloned, but passed by reference.
	 * 
	 * @param toolName
	 *            the name of the tool
	 * @param parameters
	 *            the {@link Parameter}s
	 */
	public ToolParameterSet( String toolName, Parameter... parameters) {
		super(parameters);
		this.toolName = toolName;
	}
	
	/**
	 * Constructs a {@link ToolParameterSet} given a tool name and some {@link Parameter}s.
	 * The {@link Parameter}s are not cloned, but passed by reference.
	 * 
	 * @param toolName
	 *            the name of the tool
	 * @param parameters
	 *            the {@link Parameter}s
	 */
	public ToolParameterSet( String toolName, Collection<Parameter> parameters) {
		super(parameters);
		this.toolName = toolName;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ToolParameterSet} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ToolParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>representation</code>
	 */
	public ToolParameterSet(StringBuffer representation) throws NonParsableException {
		super(representation);
	}
	
	/**
	 * Returns the name of the tool that can be run with this {@link ToolParameterSet}.
	 *  
	 * @return the name of the tool
	 */
	public String getToolName() {
		return toolName;
	}
	
	public ToolParameterSet clone() throws CloneNotSupportedException {
		return (ToolParameterSet) super.clone();
	}
	
	public StringBuffer toXML() {
		StringBuffer representation = new StringBuffer();
		XMLParser.appendObjectWithTags(representation, toolName, "toolName");
		representation.append(super.toXML());
		XMLParser.addTags(representation, "ToolParameterSet");
		return representation;
	}
	protected void fromXML( StringBuffer representation ) throws NonParsableException {
		representation = XMLParser.extractForTag(representation, "ToolParameterSet");
		toolName = (String) XMLParser.extractObjectForTags(representation, "toolName");
		super.fromXML(representation);
	}
	
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine, int indentation ) throws Exception {
		if( getParent()!=null ) {
			super.toGalaxy(namePrefix, configPrefix + toolName +".", depth, descBuffer, configBuffer, addLine, XMLParser.nextIndentation(indentation));
			XMLParser.addTagsAndAttributes(descBuffer, "section", "name=\"" + toolName + "\" title=\"" + toolName + " parameters\" expanded=\"" + !hasDefaultOrIsSet() + "\"", indentation );
			descBuffer.insert(0, "\n");
			descBuffer.append("\n");
		} else {
			super.toGalaxy(namePrefix, configPrefix, depth, descBuffer, configBuffer, addLine, indentation);
		}
	}
	
	public void toGalaxyTest( String namePrefix, int depth, StringBuffer testBuffer, int indentation ) throws Exception {
		if( getParent()!=null ) {
			StringBuffer help = new StringBuffer();
			super.toGalaxyTest(namePrefix, depth, help, XMLParser.nextIndentation(indentation));
			if( help.length()>0 ) {
				XMLParser.addTagsAndAttributes(help, "section", "name=\"" + toolName + "\"", indentation );
				help.insert(0, "\n");
				help.append("\n");
				testBuffer.append(help);
			}
		} else {
			super.toGalaxyTest(namePrefix, depth, testBuffer, indentation);
		}
	}
}
