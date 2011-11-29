/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.parameters;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.utils.galaxy.GalaxyAdaptor;

/**
 * Class for a {@link ParameterSetContainer} that contains a
 * {@link ParameterSet} as value. This {@link ParameterSet} can be set either
 * using the constructor or using the method {@link #setValue(Object)} and can
 * be obtained by the method {@link #getValue()}. {@link ParameterSetContainer}s
 * can be used to build tree-structures of {@link Parameter}s.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ParameterSetContainer extends Parameter implements Rangeable,
		RangeIterator, GalaxyConvertible {

	private String name;
	private String comment;
	private ParameterSet parameters;

	/**
	 * The message of the last error or <code>null</code> if no error occurred.
	 */
	protected String errorMessage;

	/**
	 * Creates an new {@link ParameterSetContainer} out of a
	 * {@link ParameterSet}.
	 * 
	 * @param name
	 *            the name of the {@link ParameterSetContainer}
	 * @param comment
	 *            a comment on the {@link ParameterSetContainer}
	 * @param content
	 *            the content of the {@link ParameterSetContainer} (the
	 *            contained {@link ParameterSet})
	 */
	public ParameterSetContainer(String name, String comment,
			ParameterSet content) {
		this.name = name;
		this.comment = comment;
		this.parameters = content;
		this.parameters.setParent(this);
		this.errorMessage = null;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ParameterSetContainer} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public ParameterSetContainer(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public ParameterSetContainer clone() throws CloneNotSupportedException {
		ParameterSetContainer clone = (ParameterSetContainer) super.clone();
		clone.parameters = parameters == null ? null : parameters.clone();
		if (clone.parameters != null) {
			clone.parameters.setParent(clone);
		}
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getName()
	 */
	@Override
	public String getName() {
		return name;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getDatatype()
	 */
	@Override
	public DataType getDatatype() {
		return DataType.PARAMETERSET;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getComment()
	 */
	@Override
	public String getComment() {
		return comment;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isRequired()
	 */
	@Override
	public boolean isRequired() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#checkValue(java.lang.Object)
	 */
	@Override
	public boolean checkValue(Object value) {
		return (value == null || this.parameters.getClass().isInstance( value  ) );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setValue(java.lang.Object)
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		if (checkValue(value)) {
			this.parameters = (ParameterSet) value;
			this.parameters.setParent(this);
		} else {
			throw new IllegalValueException(
					"Only parameter sets allowed for ParameterSetContainer!");
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getValue()
	 */
	@Override
	public ParameterSet getValue() {
		return parameters;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		return parameters != null && parameters.hasDefaultOrIsSet();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isSet()
	 */
	@Override
	public boolean isSet() {
		return parameters != null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isAtomic()
	 */
	@Override
	public boolean isAtomic() {
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getErrorMessage()
	 */
	@Override
	public String getErrorMessage() {
		return parameters.getErrorMessage();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#simplify()
	 */
	@Override
	public void simplify() {
		parameters.simplify();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#reset()
	 */
	@Override
	public void reset() {
		parameters.reset();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws Exception {
		throw new Exception("Not applicable to ParameterSetContainer");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags(buf, "superParameter");
		XMLParser.appendObjectWithTags(buf, name, "name");
		XMLParser.appendObjectWithTags(buf, comment, "comment");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, parameters, "parameters");
		XMLParser.addTags(buf, "parameterSetContainer");
		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,"parameterSetContainer");
		super.fromXML(XMLParser.extractForTag(representation,"superParameter"));
		name = XMLParser.extractObjectForTags(representation, "name", String.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		comment = XMLParser.extractObjectForTags(representation, "comment", String.class );
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(representation, "errorMessage", String.class ) );
		parameters = XMLParser.extractObjectForTags( representation, "parameters", ParameterSet.class );
		parameters.setParent(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#isRangeable()
	 */
	public boolean isRangeable() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#getRangedInstance()
	 */
	public Parameter getRangedInstance() throws Exception {
		this.parameters.makeRanged();
		return this;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#next()
	 */
	public boolean next() throws ParameterException {
		return parameters.next();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#resetToFirst()
	 */
	public void resetToFirst() {
		parameters.resetToFirst();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#isRanged()
	 */
	public boolean isRanged() {
		return parameters.isRanged();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#getNumberOfValues()
	 */
	public int getNumberOfValues() {
		return parameters.getNumberOfValues();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#valuesToString()
	 */
	public String valuesToString() {
		// return name;
		return "(" + parameters.valuesToString() + ")";
	}

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer ) throws Exception {
		StringBuffer pars = new StringBuffer();
		((GalaxyConvertible)parameters).toGalaxy( namePrefix, configPrefix, depth+1, pars, configBuffer );
		String color = GalaxyAdaptor.getColor( depth );
		
		StringBuffer buf = new StringBuffer();
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"hidden\" name=\""+namePrefix+"_contbegin"+"\" help=\"&lt;hr style=&quot;height:2px;background-color:"+color+";color:"+color+";border:none&quot; /&gt;\"" );
		buf.append( pars );
		StringBuffer buf2 = new StringBuffer();
		XMLParser.addTagsAndAttributes( buf2, "param", "type=\"hidden\" name=\""+namePrefix+"_contend"+"\" help=\"&lt;hr style=&quot;height:2px;background-color:"+color+";color:"+color+";border:none&quot; /&gt;\"" );
		buf.append( buf2 );
		
		descBuffer.append( buf );
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		((GalaxyConvertible)parameters).fromGalaxy( namePrefix, command );	
	}
	
	

}
