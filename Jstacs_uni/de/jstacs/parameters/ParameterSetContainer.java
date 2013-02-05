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
import de.jstacs.io.NonParsableException;
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
public class ParameterSetContainer extends Parameter implements GalaxyConvertible {

	private ParameterSet parameters;
	private Class<? extends ParameterSet> parameterClass;

	/**
	 * The message of the last error or <code>null</code> if no error occurred.
	 */
	protected String errorMessage;

	/**
	 * Creates an new {@link ParameterSetContainer} out of a {@link ParameterSet}.
	 * 
	 * @param p
	 *            the content of the {@link ParameterSetContainer} (the
	 *            contained {@link ParameterSet})
	 *            
	 * @see ParameterSet#getName(ParameterSet)
	 * @see ParameterSet#getComment(ParameterSet)          
	 */
	public ParameterSetContainer( ParameterSet p ) {
		this( ParameterSet.getName(p), ParameterSet.getComment(p), p );
	}
	
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
	public ParameterSetContainer(String name, String comment, ParameterSet content) {
		this( name, comment, content.getClass() );
		this.parameters = content;
		this.parameters.setParent(this);
	}
	
	/**
	 * Creates an new {@link ParameterSetContainer} out of the class
	 * of a {@link ParameterSet}.
	 * 
	 * @param contentClazz
	 *            the class of the contained {@link ParameterSet}
	 *            
	 * @see ParameterSet#getName(Class)
	 * @see ParameterSet#getComment(Class)
	 */
	public ParameterSetContainer( Class<? extends ParameterSet> contentClazz) {
		this( ParameterSet.getName(contentClazz), ParameterSet.getComment(contentClazz), contentClazz );
	}

	
	/**
	 * Creates an new {@link ParameterSetContainer} out of the class
	 * of a {@link ParameterSet}.
	 * 
	 * @param name
	 *            the name of the {@link ParameterSetContainer}
	 * @param comment
	 *            a comment on the {@link ParameterSetContainer}
	 * @param contentClazz
	 *            the class of the contained {@link ParameterSet}
	 */
	public ParameterSetContainer(String name, String comment, Class<? extends ParameterSet> contentClazz) {
		super( name, comment, DataType.PARAMETERSET );
		this.parameterClass = contentClazz;
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
	public ParameterSetContainer(StringBuffer representation) throws NonParsableException {
		super( representation );
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
		if( value == null || !(value instanceof ParameterSet) ) {
			return false;
		}
		return (this.parameterClass.isInstance( value ) || ( parameters==null || parameters.isComparable( (ParameterSet) value ) ) );
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
			System.out.println( value.getClass().getSimpleName() );
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
		if(parameters == null){
			loadParameters();
		}
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
		if(parameters == null){
			return null;
		}else{
			return parameters.getErrorMessage();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#reset()
	 */
	@Override
	public void reset() {
		if(parameters != null){
			parameters.reset();
		}
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
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "parameterSetContainer";
	}

	private void loadParameters() {
		if(parameters == null){
			try{
				parameters = parameterClass.newInstance();
				parameters.setParent( this );
			}catch(Exception e){
				RuntimeException ex = new RuntimeException( e );
				ex.setStackTrace( e.getStackTrace() );
				throw ex;
			}
			parameters.setParent(this);
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		if(parameters == null){
			XMLParser.appendObjectWithTags(buf, parameterClass, "parameterClass");
		}else{
			XMLParser.appendObjectWithTags(buf, parameters, "parameters");
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		super.appendFurtherInfos( representation );
		
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(representation, "errorMessage", String.class ) );
		if( !XMLParser.hasTag(representation, "parameters", null, null) ){
			this.parameterClass = XMLParser.extractObjectForTags( representation, "parameterClass", Class.class );
		}else{
			parameters = XMLParser.extractObjectForTags( representation, "parameters", ParameterSet.class );
			this.parameterClass = parameters.getClass();
			parameters.setParent(this);
		}
	}

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine ) throws Exception {
	
		if(parameters == null){
			loadParameters();
		}
		StringBuffer pars = new StringBuffer();
		((GalaxyConvertible)parameters).toGalaxy( namePrefix, configPrefix, depth, pars, configBuffer, false );
		//String color = GalaxyAdaptor.getColor( depth );
		
	/*	StringBuffer buf = new StringBuffer();
		//XMLParser.addTagsAndAttributes( buf, "param", "type=\"hidden\" name=\""+namePrefix+"_contbegin"+"\" help=\"&lt;hr style=&quot;height:2px;background-color:"+color+";color:"+color+";border:none&quot; /&gt;\"" );
		buf.append( pars );
		StringBuffer buf2 = new StringBuffer();
		//XMLParser.addTagsAndAttributes( buf2, "param", "type=\"hidden\" name=\""+namePrefix+"_contend"+"\" help=\"&lt;hr style=&quot;height:2px;background-color:"+color+";color:"+color+";border:none&quot; /&gt;\"" );
		buf.append( buf2 );*/
		
		descBuffer.append( pars );
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		if(parameters == null){
			loadParameters();
		}
		((GalaxyConvertible)parameters).fromGalaxy( namePrefix, command );	
	}
	
	

}
