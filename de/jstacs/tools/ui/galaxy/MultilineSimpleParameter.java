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

package de.jstacs.tools.ui.galaxy;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.GalaxyConvertible;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.ParameterValidator;

/**
 * An extension of {@link SimpleParameter} that renders as a textarea in Galaxy, which is only suitable for {@link DataType#STRING}s.
 * Besides the {@link MultilineSimpleParameter#toGalaxy( String, String, int, StringBuffer, StringBuffer, boolean, int)}
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
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine, int indentation ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		StringBuffer buf = new StringBuffer();
		if(validator != null && validator instanceof GalaxyConvertible){
			((GalaxyConvertible)validator).toGalaxy( namePrefix+"_valid", null, depth, buf, null, false, XMLParser.nextIndentation(indentation) );
		}
		String line = "";
		if(addLine){
			line = "&lt;hr /&gt;";
		}
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"text\" area=\"true\" size=\"10x80\" name=\""+namePrefix+"\" label=\""+line+getName()+"\" help=\""+getComment()+"\" value=\""+(defaultValue == null ? "" : defaultValue)+"\" optional=\""+(!isRequired())+"\"", indentation );
		descBuffer.append( buf );
		
		buf = new StringBuffer();
		buf.append( "${"+configPrefix+namePrefix+"}" );
		XMLParser.addTags( buf, namePrefix );
		configBuffer.append( buf );		
	}	
}
