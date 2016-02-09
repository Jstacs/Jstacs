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

import java.util.HashSet;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;

/**
 * Class for a collection parameter, i.e. a parameter that provides some
 * collection of possible values the user can choose from. Instances of this
 * class can be used as parameters in a {@link ParameterSet}.
 * 
 * @see de.jstacs.parameters.ParameterSet
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public abstract class AbstractSelectionParameter extends Parameter implements Rangeable, GalaxyConvertible {

	/**
	 * The internal {@link ParameterSet} that holds the possible values
	 */
	protected ParameterSet parameters;

	/**
	 * <code>true</code> if the user has selected an item
	 */
	protected boolean userSelected;

	/**
	 * <code>true</code> if the current parameter is required, false otherwise
	 */
	private boolean required;

	/**
	 * If a value was illegal for the collection parameter, this field holds the
	 * error message.
	 */
	protected String errorMessage;

	/**
	 * Indicates if this {@link AbstractSelectionParameter} shall be rangeable
	 */
	private boolean rangeable;

	// default constructor
	private AbstractSelectionParameter( DataType datatype, String name, String comment, boolean required ) {
		super( name, comment, datatype );
		this.required = required;
		this.userSelected = false;
		this.rangeable = true;
	}

	/**
	 * Constructor for a {@link AbstractSelectionParameter} of {@link SimpleParameter}s.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if this {@link AbstractSelectionParameter} is
	 *            required, <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the length of <code>keys</code> and the
	 *             <code>values</code> is different or the collection is
	 *             inconsistent for some other reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 * 
	 * @see #AbstractSelectionParameter(DataType, String[], Object[], String[], String, String, boolean)
	 * @see SimpleParameter
	 */
	public AbstractSelectionParameter(DataType datatype, String[] keys, Object[] values, String name, String comment, boolean required)
			throws InconsistentCollectionException, IllegalValueException,
			DatatypeNotValidException {
		this( datatype, keys, values, null, name, comment, required );
	}

	/**
	 * Constructor for a {@link AbstractSelectionParameter}.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param comments
	 *            the comments on the possible values
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the {@link AbstractSelectionParameter} is
	 *            required, <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the length of <code>keys</code> and <code>values</code> is
	 *             different or the collection is inconsistent for some other
	 *             reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 * 
	 * @see #createParameterSet(Object[], String[], String[])
	 */
	public AbstractSelectionParameter(DataType datatype, String[] keys, Object[] values, String[] comments, String name, String comment,
			boolean required) throws InconsistentCollectionException, IllegalValueException, DatatypeNotValidException {
		this(datatype, name, comment, required);

		if ( !(values instanceof Parameter[]) && (keys == null || keys.length != values.length || (comments != null && keys.length != comments.length)) ) {
			throw new InconsistentCollectionException( "You have to define the same number of keys and values for a AbstractSelectionParameter!");
		}

		createParameterSet(values, keys, comments);
	}

	/**
	 * Constructor for a {@link AbstractSelectionParameter} from an array of
	 * {@link ParameterSet}s. This constructor can be used to easily construct a
	 * {@link AbstractSelectionParameter} that lets the user select from a list of
	 * possible options that all require an own set of {@link Parameter}s to be
	 * instantiated.
	 * 
	 * @param values
	 *            the array of {@link ParameterSet}s
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the parameter is required,
	 *            <code>false</code> otherwise
	 * 
	 * @see #createParameterSet(Object[], String[], String[])
	 * @see ParameterSet#getName(ParameterSet)
	 * @see ParameterSet#getComment(ParameterSet)
	 */
	public AbstractSelectionParameter( String name, String comment, boolean required, ParameterSet... values) {
		this(DataType.PARAMETERSET, name, comment, required);
		try {
			createParameterSet(values, null, null);
		} catch( Exception shouldNotHappen ) {
			throw new RuntimeException( shouldNotHappen );
		}
	}

	/**
	 * Constructor for a {@link AbstractSelectionParameter} from an array of
	 * {@link Class}es of {@link ParameterSet}s. This constructor can be used to easily construct a
	 * {@link AbstractSelectionParameter} that lets the user select from a list of
	 * possible options that all require an own set of {@link Parameter}s to be
	 * instantiated. It is the lazy evaluation variant of {@link #AbstractSelectionParameter(String, String, boolean, ParameterSet...)},
	 * i.e., the {@link ParameterSet}s are only created if necessary.
	 * 
	 * @param values
	 *            the array of {@link Class}es of {@link ParameterSet}s
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the parameter is required,
	 *            <code>false</code> otherwise
	 * 
	 * @see #createParameterSet(Object[], String[], String[])
	 * @see ParameterSet#getName(ParameterSet)
	 * @see ParameterSet#getComment(ParameterSet)
	 */
	public AbstractSelectionParameter( String name, String comment, boolean required, Class<? extends ParameterSet>... values) {
		this(DataType.PARAMETERSET, name, comment, required);
		try {
			createParameterSet(values, null, null);
		} catch( Exception shouldNotHappen ) {
			throw new RuntimeException( shouldNotHappen );
		}
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores an instance of {@link AbstractSelectionParameter} from a XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public AbstractSelectionParameter(StringBuffer representation) throws NonParsableException {
		super(representation);
	}

	/**
	 * Creates a new {@link ParameterSet} from an array of values, an array of
	 * names and an array of comments.
	 * 
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param comments
	 *            the comments on the possible values
	 * 
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the data type is not allowed, i.e. is not one of the
	 *             primitive data types, {@link String} or
	 *             {@link DataType#PARAMETERSET}
	 * @throws InconsistentCollectionException
	 *             if the internally used keys are not unique 
	 */
	@SuppressWarnings("unchecked")
	protected void createParameterSet(Object[] values, String[] keys, String[] comments) throws DatatypeNotValidException, IllegalValueException, InconsistentCollectionException {
		Parameter[] pars = new Parameter[values.length];
		HashSet<String> hash = new HashSet<String>();
		for (int i = 0; i < pars.length; i++) {
			if( values[i] instanceof Parameter ) {
				pars[i] = (Parameter) values[i];
			} else if( values[i] instanceof ParameterSet ) {
				if( keys == null || keys[i] == null ) {
					pars[i] = new ParameterSetContainer( (ParameterSet) values[i] );
				} else {
					pars[i] = new ParameterSetContainer( keys[i], comments==null ? null : comments[i], (ParameterSet) values[i] );
				}
			} else if( values[i] instanceof Class && ParameterSet.class.isAssignableFrom( (Class) values[i] )  ) {
				if( keys == null || keys[i] == null ) {
					pars[i] = new ParameterSetContainer( (Class<? extends ParameterSet>) values[i] );
				} else {
					pars[i] = new ParameterSetContainer( keys[i], comments==null ? null : comments[i], (Class<? extends ParameterSet>) values[i] );
				}
			} else {
				if (keys == null || keys[i] == null) {
					throw new IllegalArgumentException( "You have to state the key for entity " + i);
				}
				pars[i] = new SimpleParameter(datatype, keys[i], comments==null ? null : comments[i] , false);
				pars[i].setValue(values[i]);
			}
			if( !hash.contains( pars[i].getName() ) ) {
				hash.add( pars[i].getName() );
			} else {
				throw new InconsistentCollectionException( "The key \"" + pars[i].getName() +"\" is used multiple times." );
			}
		}

		parameters = new SimpleParameterSet(pars);
		parameters.setParent(this);
		try {
			setDefault( pars[0].getName() );
		} catch (Exception doesNotHappen) {
			System.out.println( pars[0] );
			throw new RuntimeException( doesNotHappen );
		}
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public AbstractSelectionParameter clone() throws CloneNotSupportedException {
		AbstractSelectionParameter clone = (AbstractSelectionParameter) super.clone();
		clone.parameters = parameters == null ? null : parameters.clone();
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isAtomic()
	 */
	@Override
	public boolean isAtomic() {
		// return !(getValue() instanceof ParameterSet);
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (parameters.getParameterAt(i).getValue() instanceof ParameterSet) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns <code>true</code>, if this {@link AbstractSelectionParameter} has a
	 * default value.
	 * 
	 * @return if this {@link AbstractSelectionParameter} has a default value
	 * 
	 * @see #setDefault(Object)
	 */
	public abstract boolean hasDefault();

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isSet()
	 */
	@Override
	public boolean isSet() {
		return isUserSelected();
	}



	/**
	 * Returns the possible values in this collection.
	 * 
	 * @return the values/parameters in this collection
	 */
	public ParameterSet getParametersInCollection() {
		return parameters;
	}

	/**
	 * Sets the value returned by {@link #isRangeable()} to
	 * <code>rangeable</code>.
	 * 
	 * @param rangeable
	 *            the new value
	 */
	public void setRangeable(boolean rangeable) {
		this.rangeable = rangeable;
	}

	
	public boolean isComparable( Parameter p ) {
		boolean res = getClass().equals(p.getClass()) && getDatatype() == p.getDatatype() && getName().equals(p.getName()) && getComment().equals( p.getComment() );
		if(res){
			AbstractSelectionParameter sp2 = (AbstractSelectionParameter)p;
			res &= this.parameters.getNumberOfParameters() == sp2.parameters.getNumberOfParameters();
			if(res){
				for(int i=0;i<parameters.getNumberOfParameters();i++){
					res &= parameters.getParameterAt( i ).isComparable( sp2.parameters.getParameterAt( i ) );
				}
			}
		}
		return res;
	}
	
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#isRangeable()
	 */
	public boolean isRangeable() {
		return rangeable;
	}
	
	/**
	 * MSPDTest whether a given <code>value</code> can be used in {@link #setValue(Object)}.
	 * If so, the index of corresponding option is returned, otherwise -1.
	 * 
	 * @param value the value to be checked
	 * 
	 * @return the index of the corresponding option or -1 if there's no corresponding option
	 * 
	 * @see #errorMessage
	 * @see #checkValue(Object)
	 */
	@SuppressWarnings("unchecked")
	protected int check(Object value) {
		if (value == null) {
			return -1;
		}
		
		Object val2 = value;
		if( value instanceof ParameterSet ) {
			val2 = ParameterSet.getName( (ParameterSet)value );
		} else if( value instanceof Class && ParameterSet.class.isAssignableFrom( (Class) value )  ) {
			val2 = ParameterSet.getName( (Class<? extends ParameterSet>)value );
		} else if( value instanceof Parameter ) {
			val2 = ((Parameter)value).getName();
		}
		
		if ( val2 instanceof String ) {
			for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
				if (parameters.getParameterAt(i).getName().equals(val2)) {
					errorMessage = null;
					return i;
				}
			}
		}
		errorMessage = "The value is not in the set of defined values: " + value + ".";
		return -1;
	}

	/**
	 * Returns <code>true</code> if the key specified by <code>value</code> is
	 * in the set of keys of this {@link AbstractSelectionParameter}.
	 * 
	 * @return <code>true</code> if the key specified by <code>value</code> is
	 *         in the set of keys of this {@link AbstractSelectionParameter},
	 *         <code>false</code> otherwise
	 */
	@Override
	public boolean checkValue(Object value) {
		return check(value) >= 0 ;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getErrorMessage()
	 */
	@Override
	public String getErrorMessage() {
		return errorMessage;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, required, "required");
		XMLParser.appendObjectWithTags(buf, userSelected, "userSelected");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, rangeable, "rangeable");
		XMLParser.appendObjectWithTags(buf, parameters, "collection");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		super.extractFurtherInfos( representation );
		
		required = XMLParser.extractObjectForTags(representation, "required", boolean.class );
		userSelected = XMLParser.extractObjectForTags(representation, "userSelected", boolean.class );
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(representation, "errorMessage", String.class ) );
		StringBuffer help = XMLParser.extractForTag(representation, "rangeable");
		if (help == null) {
			rangeable = false;
		} else {
			rangeable = Boolean.parseBoolean(help.toString());
		}
		parameters = new SimpleParameterSet(XMLParser.extractForTag(representation,"collection"));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isRequired()
	 */
	@Override
	public boolean isRequired() {
		return required;
	}

	/**
	 * Returns <code>true</code> if the option at position <code>idx</code> is
	 * selected.
	 * 
	 * @param idx
	 *            the position
	 * 
	 * @return <code>true</code> if the option at position <code>idx</code> is
	 *         selected, <code>false</code> otherwise
	 */
	public abstract boolean isSelected(int idx);

	/**
	 * Returns <code>true</code> if the value was selected by the user.
	 * 
	 * @return <code>true</code> if the value was selected by the user.
	 */
	public boolean isUserSelected() {
		return userSelected;
	}

	/**
	 * Sets the default value of this {@link AbstractSelectionParameter} to
	 * <code>defaultValue</code>. This method also sets the current
	 * value of this {@link AbstractSelectionParameter} to the default
	 * and resets it such that {@link AbstractSelectionParameter#isUserSelected()}
	 * returns <code>false</code>.
	 * 
	 * @param defaultValue
	 *            the default value
	 * 
	 * @throws Exception
	 *             if the default value is not an appropriate value for this
	 *             {@link Parameter}
	 */
	public abstract void setDefault(Object defaultValue) throws Exception;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o2) {
		if (o2 instanceof AbstractSelectionParameter) {
			ParameterSet parSet2 = ((AbstractSelectionParameter) o2)
					.getParametersInCollection();
			if (parSet2.getNumberOfParameters() != parameters
					.getNumberOfParameters()) {
				return false;
			} else {
				if (!(((AbstractSelectionParameter) o2).getName().equals(name) && ((AbstractSelectionParameter) o2)
						.getComment().equals(comment))) {
					return false;
				}

				for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
					if (!parameters.getParameterAt(i).equals(
							parSet2.getParameterAt(i))) {
						return false;
					}
				}
				return true;
			}
		} else {
			return false;
		}
	}

	/**
	 * This exception is thrown if the {@link AbstractSelectionParameter} is
	 * inconsistent for some reason.
	 * 
	 * @author Jan Grau
	 */
	public class InconsistentCollectionException extends ParameterException {

		private static final long serialVersionUID = -2703514434545861722L;

		/**
		 * Constructs a new {@link InconsistentCollectionException} with message
		 * <code>message</code>.
		 * 
		 * @param message
		 *            the error message
		 */
		public InconsistentCollectionException(String message) {
			super(message);
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#getRangedInstance()
	 */
	@Override
	public MultiSelectionParameter getRangedInstance() throws Exception {
		//boolean[] selected = new boolean[parameters.getNumberOfParameters()];
		//selected[0] = true;
		Parameter[] p = new Parameter[parameters.getNumberOfParameters()];
		for( int i = 0; i < p.length; i++ ) {
			p[i] = parameters.getParameterAt(i);
		}
		MultiSelectionParameter par = 
			new MultiSelectionParameter( datatype, null, p, null, getName(), getComment(), required );
		// new MultiSelectionParameter( this.parameters.clone(), selected, new boolean[parameters.getNumberOfParameters()], false, getName(), getComment(), required, datatype, errorMessage, 0 );
		return par;
	}
	

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine ) throws Exception {
		StringBuffer buf = new StringBuffer();
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		
		for(int i=0;i<parameters.getNumberOfParameters();i++){
			if(isSelected( i )){
				XMLParser.appendObjectWithTagsAndAttributes( buf, parameters.getParameterAt( i ).getName() , "option", "value=\""+parameters.getParameterAt( i ).getName()+"\" selected=\"true\"", false );
			}else{
				XMLParser.appendObjectWithTagsAndAttributes( buf, parameters.getParameterAt( i ).getName() , "option", "value=\""+parameters.getParameterAt( i ).getName()+"\"", false );
			}
		}
		String line = "";
		if(addLine){
			//line = "&lt;hr style=&quot;height:2px;background-color:"+GalaxyAdaptor.getColor( depth+2 )+";color:"+GalaxyAdaptor.getColor( depth )+";border:none&quot; /&gt;";
			line = "&lt;hr /&gt;";
		}
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"select\" format=\"text\" name=\""+namePrefix+"\" label=\""+line+getName()+"\" optional=\""+(!isRequired())+"\" help=\""+getComment()/*+(isAtomic() ? "" : "&lt;hr /&gt;")*/+"\"" );
	
		StringBuffer buf3 = new StringBuffer("${"+configPrefix+(isAtomic() ? "" : namePrefix+"_cond.")+namePrefix+"}");
		XMLParser.addTags( buf3, namePrefix );
		configBuffer.append( buf3 );
		if(!isAtomic()){
			for(int i=0;i<parameters.getNumberOfParameters();i++){
				StringBuffer temp = new StringBuffer();
				StringBuffer temp2 = new StringBuffer();
				((GalaxyConvertible)parameters.getParameterAt( i )).toGalaxy( namePrefix+"_opt"+i, configPrefix+namePrefix+"_cond.", depth+1, temp, temp2, false );
				XMLParser.addTagsAndAttributes( temp, "when", "value=\""+parameters.getParameterAt( i ).getName()+"\"" );
				buf.append( temp );
				configBuffer.append( (i== 0 ? "#if " : "#elif ")+"$"+configPrefix+namePrefix+"_cond."+namePrefix+" == \""+parameters.getParameterAt( i ).getName()+"\"\n" );
				configBuffer.append( temp2 );
				configBuffer.append( "\n" );
			}
			configBuffer.append( "#end if\n" );
		}
		
		
		if(!isAtomic()){
			XMLParser.addTagsAndAttributes( buf, "conditional", "name=\""+namePrefix+"_cond\"" );
		}
		
		descBuffer.append( buf );
	}
}