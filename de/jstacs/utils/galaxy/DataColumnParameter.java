package de.jstacs.utils.galaxy;

import java.io.PrintWriter;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.GalaxyConvertible;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.ParameterValidator;


public class DataColumnParameter extends SimpleParameter {

	private String dataRef;
	
	/**
	 * @param datatype
	 * @param name
	 * @param comment
	 * @param required
	 * @param defaultVal
	 * @throws DatatypeNotValidException
	 * @throws IllegalValueException
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required, Integer defaultVal )
																														throws DatatypeNotValidException,
																														IllegalValueException {
		super( DataType.INT, name, comment, required, defaultVal );
		this.dataRef = dataRef;
	}

	/**
	 * @param datatype
	 * @param name
	 * @param comment
	 * @param required
	 * @param validator
	 * @param defaultVal
	 * @throws ParameterException
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required, ParameterValidator validator,
								Integer defaultVal ) throws ParameterException {
		super( DataType.INT, name, comment, required, validator, defaultVal );
		this.dataRef = dataRef;
	}

	/**
	 * @param datatype
	 * @param name
	 * @param comment
	 * @param required
	 * @param validator
	 * @throws DatatypeNotValidException
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required, ParameterValidator validator )
																																throws DatatypeNotValidException {
		super( DataType.INT, name, comment, required, validator );
		this.dataRef = dataRef;
	}

	/**
	 * @param datatype
	 * @param name
	 * @param comment
	 * @param required
	 * @throws DatatypeNotValidException
	 */
	public DataColumnParameter( String dataRef, String name, String comment, boolean required ) throws DatatypeNotValidException {
		super( DataType.INT, name, comment, required );
		this.dataRef = dataRef;
	}

	/**
	 * @param representation
	 * @throws NonParsableException
	 */
	public DataColumnParameter( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	protected String dataTypeToGalaxy() {
		
		return "data_column";
		
	}
	
	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine ) throws Exception {
		String refPrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( dataRef );
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
	
	
	
	
}
