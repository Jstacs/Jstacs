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

import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.AbstractSelectionParameter;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.GalaxyConvertible;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.ParameterValidator;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.utils.Pair;

/**
 * {@link SimpleParameter} that represents a data column parameter in Galaxy and JstacsFX.
 * 
 * @author Jan Grau
 *
 */
public class DataColumnParameter extends SimpleParameter {

	private String dataRef;
	
	/**
	 * Creates a new {@link DataColumnParameter} with given name, comment, and reference.
	 * @param dataRef the ID of the referenced parameter (tabular) in Galaxy.
	 * @param name the name of the parameter
	 * @param comment a comment on the parameter
	 * @param required if this parameter is required
	 * @param defaultVal the default value
	 * @throws DatatypeNotValidException should not happen
	 * @throws IllegalValueException should not happen
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required, Integer defaultVal ) throws DatatypeNotValidException,
																														IllegalValueException {
		super( DataType.INT, name, comment, required, defaultVal );
		this.dataRef = dataRef;
	}

	/**
	 * Creates a new {@link DataColumnParameter} with given name, comment, and reference.
	 * @param dataRef the ID of the referenced parameter (tabular) in Galaxy.
	 * @param name the name of the parameter
	 * @param comment a comment on the parameter
	 * @param required if this parameter is required
	 * @param validator a validator of admissible values
	 * @param defaultVal the default value
	 * @throws ParameterException if the default value is not admissible
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required, ParameterValidator validator,
								Integer defaultVal ) throws ParameterException {
		super( DataType.INT, name, comment, required, validator, defaultVal );
		this.dataRef = dataRef;
	}

	/**
	 * Creates a new {@link DataColumnParameter} with given name, comment, and reference.
	 * @param dataRef the ID of the referenced parameter (tabular) in Galaxy.
	 * @param name the name of the parameter
	 * @param comment a comment on the parameter
	 * @param required if this parameter is required
	 * @param validator a validator of admissible values
	 * @throws DatatypeNotValidException should not happen
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required, ParameterValidator validator ) throws DatatypeNotValidException {
		super( DataType.INT, name, comment, required, validator );
		this.dataRef = dataRef;
	}

	/**
	 * Creates a new {@link DataColumnParameter} with given name, comment, and reference.
	 * @param dataRef the ID of the referenced parameter (tabular) in Galaxy.
	 * @param name the name of the parameter
	 * @param comment a comment on the parameter
	 * @param required if this parameter is required
	 * @throws DatatypeNotValidException should not happen
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required ) throws DatatypeNotValidException {
		super( DataType.INT, name, comment, required );
		this.dataRef = dataRef;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link DataColumnParameter} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DataColumnParameter} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public DataColumnParameter( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	protected String dataTypeToGalaxy() {
		
		return "data_column";
		
	}
	
	public String getDataRef(){
		return dataRef;
	}
	
	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine ) throws Exception {
		//String refPrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( dataRef );
		
		FileParameter fp = find(this.getParent(), dataRef);
		
		String tool = namePrefix.substring(0,namePrefix.indexOf("_ps"));
		
		String refPrefix = tool+getNamePostfix(fp);
		
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		StringBuffer buf = new StringBuffer();
		if(validator != null && validator instanceof GalaxyConvertible){
			((GalaxyConvertible)validator).toGalaxy( namePrefix+"_valid", null, depth, buf, null, false );
		}
		
		String line = "";
		if(addLine){
			line = "&lt;hr /&gt;";
		}
		
		XMLParser.addTagsAndAttributes( buf, "param", "type=\""+dataTypeToGalaxy()+"\" name=\""+namePrefix+"\" data_ref=\""+refPrefix+"\" force_select=\""+isRequired()+"\"  label=\""+line+getName()+"\" help=\""+getComment()+"\" value=\""+(defaultValue == null ? "" : defaultValue)+"\" optional=\""+(!isRequired())+"\"" );
		descBuffer.append( buf );
		
		buf = new StringBuffer();
		buf.append( "${"+configPrefix+namePrefix+"}" );
		XMLParser.addTags( buf, namePrefix );
		configBuffer.append( buf );		
	}
	
	private String getNamePostfix(FileParameter par) {
		String postfix = "_"+GalaxyAdaptor.getLegalName(par.getName());
		
		return getNamePostfix(par,postfix);
		
	}
	
	private String getNamePostfix(Parameter par, String postfix){
		ParameterSet set = par.getParent();
		Parameter parent = set.getParent();
		
		if(parent != null && parent instanceof AbstractSelectionParameter){
			
			ParameterSet in = ((AbstractSelectionParameter)parent).getParametersInCollection();
			for(int i=0;i<in.getNumberOfParameters();i++){
				ParameterSetContainer cont = (ParameterSetContainer) in.getParameterAt(i);
				if(cont.getValue() == set){
					return getNamePostfix(parent,"_"+GalaxyAdaptor.getLegalName(parent.getName())+"_opt"+i+postfix);
				}
			}
			
		}else{
			postfix = "_ps"+postfix;
			if(parent != null){
				postfix = getNamePostfix(parent,postfix);
			}
			return postfix;
		}
		return null;
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		
		try{
			String val = XMLParser.extractForTag( command, namePrefix ).toString();
			if(!"None".equals( val )){
				this.setValue( val );
			}
		}catch(NullPointerException e){
			throw new NullPointerException( getName()+" "+command+" "+namePrefix );
		}
	}
	
	
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, dataRef, "dataRef");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		super.extractFurtherInfos( representation );
		dataRef = (String)XMLParser.extractObjectForTags( representation, "dataRef" );
	}
	
	
	
	
	public static FileParameter find(ParameterSet top, String dataRef) {
		
		while(top.getParent() != null){
			top = top.getParent().getParent();
		}
		
		LinkedList<ParameterSet> next = new LinkedList<ParameterSet>();
		for(int i=0;i<top.getNumberOfParameters();i++){
			Parameter par = top.getParameterAt(i);
			//System.out.println(par.getName()+" "+dataRef+" "+par.getClass()+" "+(par instanceof FileParameter)+" "+(par.getName().equals(dataRef)));
			if(par instanceof FileParameter && par.getName().equals(dataRef)){
				//System.out.println("return");
				return (FileParameter)par;
			}
			if(par instanceof ParameterSetContainer){
				next.add(((ParameterSetContainer)par).getValue());
			}
		}
		for(int i=0;i<next.size();i++){
			FileParameter par = find(next.get(i),dataRef);
			if(par != null){
				return par;
			}
		}
		return null;
	}
	
	
	
	
}
