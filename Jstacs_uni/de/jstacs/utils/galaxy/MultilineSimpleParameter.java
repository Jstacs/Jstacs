package de.jstacs.utils.galaxy;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.GalaxyConvertible;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.ParameterValidator;

/**
 * An extension of {@link SimpleParameter} that renders as a textarea in Galaxy, which is only suitable for {@link DataType#STRING}s.
 * Besides the {@link MultilineSimpleParameter#toGalaxy(String, String, int, StringBuffer, StringBuffer, boolean)}
 * and {@link MultilineSimpleParameter#fromGalaxy(String, StringBuffer)}, all functionality is inherited from {@link SimpleParameter}.
 * 
 * @author Jan Grau
 */
public class MultilineSimpleParameter extends SimpleParameter {

	/**
	 * Creates a new {@link MultilineSimpleParameter} with given default value.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * @param defaultVal
	 *            the default value
	 * 
	 * @throws DatatypeNotValidException if <code>datatype</code> is not in
	 *             the values allowed for a {@link SimpleParameter}
	 * @throws IllegalValueException
	 *             if the default value is not a valid value with respect
	 *             to <code>datatype</code>
	 */
	public MultilineSimpleParameter( String name, String comment, boolean required, Object defaultVal )throws DatatypeNotValidException, IllegalValueException {
		super( DataType.STRING, name, comment, required, defaultVal );
	}

	/**
	 * Creates a new {@link MultilineSimpleParameter} with given default value and a {@link ParameterValidator}.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * @param validator
	 *            the validator for the parameter values
	 * @param defaultVal
	 *            the default value
	 * 
	 * @throws ParameterException
	 *             if either the default value is not a valid value with respect
	 *             to <code>datatype</code> or <code>datatype</code> is not in
	 *             the values allowed for a {@link SimpleParameter}
	 */
	public MultilineSimpleParameter( String name, String comment, boolean required, ParameterValidator validator,
										Object defaultVal ) throws ParameterException {
		super( DataType.STRING, name, comment, required, validator, defaultVal );
		// TODO Auto-generated constructor stub
	}

	/**
	 * Creates a new {@link MultilineSimpleParameter} with no default value and a {@link ParameterValidator}.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * @param validator
	 *            the validator for the parameter values
	 * 
	 * @throws DatatypeNotValidException
	 *             if <code>datatype</code> is not in the values allowed for a
	 *             {@link SimpleParameter}
	 */
	public MultilineSimpleParameter( String name, String comment, boolean required, ParameterValidator validator ) throws DatatypeNotValidException {
		super( DataType.STRING, name, comment, required, validator );
		// TODO Auto-generated constructor stub
	}

	/**
	 * Creates a new {@link MultilineSimpleParameter} with no default value.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * 
	 * @throws DatatypeNotValidException
	 *             if <code>datatype</code> is not one of the allowed
	 *             {@link DataType}s
	 */
	public MultilineSimpleParameter( String name, String comment, boolean required ) throws DatatypeNotValidException {
		super( DataType.STRING, name, comment, required );
		// TODO Auto-generated constructor stub
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MultilineSimpleParameter} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SimpleParameter} could not be restored from the
	 *             {@link StringBuffer} <code>representation</code>
	 */
	public MultilineSimpleParameter( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		StringBuffer buf = new StringBuffer();
		if(validator != null && validator instanceof GalaxyConvertible){
			((GalaxyConvertible)validator).toGalaxy( namePrefix+"_valid", null, depth, buf, null, false );
		}
		String line = "";
		if(addLine){
			line = "&lt;hr /&gt;";
		}
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"text\" area=\"true\" size=\"10x80\" name=\""+namePrefix+"\" label=\""+line+getName()+"\" help=\""+getComment()+"\" value=\""+(defaultValue == null ? "" : defaultValue)+"\" optional=\""+(!isRequired())+"\"" );
		descBuffer.append( buf );
		
		buf = new StringBuffer();
		buf.append( "${"+configPrefix+namePrefix+"}" );
		XMLParser.addTags( buf, namePrefix );
		configBuffer.append( buf );		
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		String val = XMLParser.extractForTag( command, namePrefix ).toString();
		
		val = unescape( val );
		this.setValue( val );
	}

	
	private static String[][] table = {
	                                   { "&", "__amp__" },
	                                   { "\"", "__quot__" },
	                                   { "'", "__apos__" },
	                                   { ">", "__gt__" },
	                                   { "<", "__lt__" },
	                                   { "\n", "(__cn__|__cr__)+"},
	                                   { "[", "__ob__" },
	                                   { "]", "__cb__" }
	};



	private static String unescape( String original ) {
		//System.out.println("before: "+original);
		if( original != null ) {
			for( int i = 0; i < table.length; i++ ) {
				original = original.replaceAll( table[i][1], table[i][0] );
			}
		}
		//System.out.println("after: "+original);
		return original;
	}

	
}
