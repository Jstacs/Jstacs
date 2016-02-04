/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.tools;

import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.AbstractSelectionParameter;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.ListResult;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;

/**
 * Class for the results of a {@link JstacsTool}.
 * In addition to a standard {@link ResultSetResult}, this class stores some supplementary information
 * like the parameters of the tool's run, the name of the tool, and the time the tool's run finished.
 * 
 * @author Jan Grau
 *
 */
public class ToolResult extends ResultSetResult {

	private ParameterSet toolParameters;
	private String toolName;
	private Date finished;
	
	/**
	 * Returns the name of the tool (see {@link JstacsTool#getToolName()}) used to create these results.
	 * @return the name
	 */
	public String getToolName() {
		return toolName;
	}

	/**
	 * Creates a new {@link ToolResult} with most arguments identical to those of a {@link ListResult}.
	 * In addition, it stores the name of the creating tool, the parameters for the tool's run (see {@link JstacsTool#getToolParameters()}), 
	 * and the date/time of result creation. 
	 * @param name the name of the result
	 * @param comment a comment on the meaning of the result
	 * @param annotation optional annotation, may be <code>null</code>
	 * @param result the set of all results of the tool's run
	 * @param toolParameters the (filled) parameters of the tool's run
	 * @param toolName the name of the tool
	 * @param finished the date/time the tool's run finished and results were created
	 */
	public ToolResult( String name, String comment, ResultSet annotation, ResultSet result, ParameterSet toolParameters, String toolName, Date finished ) {
		super( name, comment, annotation, result );
		this.toolParameters = toolParameters;

		this.toolName = toolName;
		this.finished = finished;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ToolResult} from the corresponding XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer}<code>representation</code> could
	 *             not be parsed
	 */
	public ToolResult( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
	/**
	 * Returns the date and time, when the tool's run resulting in this {@link ToolResult} finished.
	 * @return the date/time
	 */
	public Date getFinishedDate(){
		return finished;
	}

	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		XMLParser.appendObjectWithTags( buf, toolParameters, "toolParameters" );
		XMLParser.appendObjectWithTags( buf, toolName, "toolName" );
		XMLParser.appendObjectWithTags( buf, finished.getTime(), "date" );
	}

	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		super.extractFurtherInfos( representation );
		toolParameters = (ParameterSet)XMLParser.extractObjectForTags( representation, "toolParameters" );
		toolName = (String)XMLParser.extractObjectForTags( representation, "toolName" );
		long date = (Long)XMLParser.extractObjectForTags( representation, "date" );
		finished = new Date( date );
	}
	
	
	/**
	 * Sets the values of all parameters in <code>other</code> to those stored in the internal parameters
	 * that have been supplied upon construction. Works successfully only if both parameter sets have an identical structure, which 
	 * should be the case if both have been created from the same tool's {@link JstacsTool#getToolParameters()} method.
	 * @param other the parameters to be filled
	 * @see #ToolResult(String, String, ResultSet, ResultSet, ParameterSet, String, Date)
	 */
	public void setFromStoredParameters(ParameterSet other){
		if(toolParameters.isComparable( other )){
			try {
				setFromStoredParameters( toolParameters, other );
			} catch ( IllegalValueException e ) {
				RuntimeException re = new RuntimeException( e );
				throw re;
			}
		}else{
			throw new RuntimeException( "ParameterSets not comparable" );
		}
	}
	
	private void setFromStoredParameters(ParameterSet mine, ParameterSet other) throws IllegalValueException{
		for(int i=0;i<mine.getNumberOfParameters();i++){
			Parameter p1 = mine.getParameterAt( i );
			Parameter p2 = other.getParameterAt( i );
			if(p1.getDatatype() != DataType.PARAMETERSET){
				if(p1 instanceof FileParameter){
					if(( (FileParameter)p1 ).getFileContents() != null){
						p2.setValue( ( (FileParameter)p1 ).getFileContents() );
					}
				}else{
					if(p1.getValue() != null){
						p2.setValue( p1.getValue() );
					}
				}
			}else if(p1 instanceof AbstractSelectionParameter){
				ParameterSet incoll = ( (AbstractSelectionParameter)p1 ).getParametersInCollection();
				LinkedList<String> name = new LinkedList<String>();
				for(int j=0;j<incoll.getNumberOfParameters();j++){
					if( ( (AbstractSelectionParameter)p1 ).isSelected( j ) ){
						name.add( incoll.getParameterAt( j ).getName() );
					}
				}
				
				if(name.size() == 1){
					p2.setValue( name.get( 0 ) );
				}else{
					p2.setValue( name.toArray( new String[0] ) );
				}
				ParameterSet incoll2 = ( (AbstractSelectionParameter)p2 ).getParametersInCollection();
				
				setFromStoredParameters( incoll, incoll2 );
			}else{
				ParameterSet ps1 = (ParameterSet)p1.getValue();
				ParameterSet ps2 = (ParameterSet)p2.getValue();
				setFromStoredParameters( ps1, ps2 );
			}
		}
	}

	/**
	 * Returns the tool's parameters that have been used to create the results stored in this {@link ToolResult}.
	 * @return the parameters
	 * @see JstacsTool#getToolParameters()
	 */
	public ParameterSet getToolParameters() {
		return toolParameters;
	}
	
	

}
