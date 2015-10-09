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
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;


public class ToolResult extends ResultSetResult {

	private ParameterSet toolParameters;
	private String toolName;
	private Date finished;
	
	public String getToolName() {
		return toolName;
	}

	public ToolResult( String name, String comment, ResultSet annotation, ResultSet result, ParameterSet toolParameters, String toolName, Date finished ) throws CloneNotSupportedException {
		super( name, comment, annotation, result );
		this.toolParameters = toolParameters;
		/*System.out.println("Created result for parameters");
		for(int i=0;i<toolParameters.getNumberOfParameters();i++){
			System.out.println(toolParameters.getParameterAt( i ));
		}*/
		this.toolName = toolName;
		this.finished = finished;
	}

	public ToolResult( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
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

	public ParameterSet getToolParameters() {
		return toolParameters;
	}
	
	

}
