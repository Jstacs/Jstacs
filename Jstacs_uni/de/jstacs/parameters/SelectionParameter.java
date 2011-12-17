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
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.utils.galaxy.GalaxyAdaptor;

/**
 * Class for a collection parameter, i.e. a parameter that provides some
 * collection of possible values the user can choose from. Instances of this
 * class can be used as parameters in a {@link ParameterSet}.
 * 
 * @see de.jstacs.parameters.ParameterSet
 * 
 * @author Jan Grau
 */
public class SelectionParameter extends AbstractCollectionParameter implements Rangeable, GalaxyConvertible {

	/**
	 * The number of the currently selected value in <code>parameters</code>
	 */
	private int selected;

	/**
	 * The number of the option selected by default
	 */
	private int defaultSelected;

	protected void init() {
		this.defaultSelected = -1;
		this.selected = 0;
	}

	/**
	 * Constructor for a {@link SelectionParameter}.
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
	 *            <code>true</code> if this {@link SelectionParameter} is
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
	 * @see AbstractCollectionParameter#AbstractCollectionParameter(DataType, String[], Object[], String[], String, String, boolean)
	 */
	public SelectionParameter(DataType datatype, String[] keys, Object[] values, String name, String comment, boolean required)
			throws AbstractCollectionParameter.InconsistentCollectionException, IllegalValueException, DatatypeNotValidException {
		this(datatype, keys, values, null, name, comment, required);
	}

	/**
	 * Constructor for a {@link SelectionParameter}.
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
	 *            <code>true</code> if the {@link SelectionParameter} is
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
	 * @see AbstractCollectionParameter#createParameterSet(Object[], String[], String[])
	 */
	public SelectionParameter(DataType datatype, String[] keys, Object[] values, String[] comments, String name, String comment,
			boolean required) throws AbstractCollectionParameter.InconsistentCollectionException, IllegalValueException, DatatypeNotValidException {
		super( datatype, keys, values, comments, name, comment, required );
	}

	/**
	 * Constructor for a {@link SelectionParameter} from an array of
	 * {@link ParameterSet}s. This constructor can be used to easily construct a
	 * {@link SelectionParameter} that lets the user select from a list of
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
	 * @see AbstractCollectionParameter#createParameterSet(Object[], String[], String[])
	 */
	public SelectionParameter( String name, String comment, boolean required, ParameterSet... values) throws DatatypeNotValidException, IllegalValueException, InconsistentCollectionException {
		super( name, comment, required, values );
	}
	
	public SelectionParameter( String name, String comment, boolean required, Class<? extends ParameterSet>... values) throws DatatypeNotValidException, IllegalValueException, InconsistentCollectionException {
		super( name, comment, required, values );
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores an instance of {@link SelectionParameter} from a XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public SelectionParameter(StringBuffer representation) throws NonParsableException {
		super(representation);
	}

	/**
	 * Returns <code>true</code>, if this {@link SelectionParameter} has a
	 * default value.
	 * 
	 * @return if this {@link SelectionParameter} has a default value
	 * 
	 * @see AbstractCollectionParameter#setDefault(Object)
	 */
	public boolean hasDefault() {
		return defaultSelected > -1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		if (getValue() instanceof ParameterSet) {
			if (!((ParameterSet) getValue()).hasDefaultOrIsSet()) {
				if (((ParameterSet) getValue()).getErrorMessage() != null) {
					errorMessage = "Selected value has the following error: "
							+ ((ParameterSet) getValue()).getErrorMessage();
				}
				return false;
			} else {
				return true;
			}
		} else if (isSet() || selected == defaultSelected) {
			return true;
		} else {
			return false;
		}
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
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "collectionParameter";
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, selected, "selected");
		XMLParser.appendObjectWithTags(buf, defaultSelected, "defaultSelected");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		super.extractFurtherInfos( representation );
		
		selected = XMLParser.extractObjectForTags(representation, "selected", int.class );
		defaultSelected = XMLParser.extractObjectForTags(representation, "defaultSelected", int.class );
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
	public boolean isSelected(int idx) {
		return idx == selected;
	}
	
	/**
	 * Returns the index of the selected value.
	 * 
	 * @return the index of the selected value
	 */
	public int getSelected() {
		return selected;
	}

	/**
	 * Returns the index of the default selected value.
	 * 
	 * @return the index of the default selected value
	 */
	protected int getDefault() {
		return defaultSelected;
	}

	/**
	 * Sets the selected value to the one that is specified by the key
	 * <code>value</code>.
	 * 
	 * @param value
	 *            the key of the desired value
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		int i = check(value);
		if( i < 0 ) {
			throw new IllegalValueException( errorMessage );
		} else {
			selected = i;
			userSelected = true;
			//TODO okay?
			if( value instanceof ParameterSetContainer ) {
				parameters.getParameterAt( selected ).setValue( ((ParameterSetContainer)value).getValue() );
			} else if( value instanceof ParameterSet ) {
				parameters.getParameterAt( selected ).setValue( value );
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws IllegalValueException {
		setValue(defaultValue);
		defaultSelected = selected;
		userSelected = false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#reset()
	 */
	@Override
	public void reset() {
		selected = defaultSelected;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getValue()
	 */
	@Override
	public Object getValue() {
		if (selected < parameters.getNumberOfParameters()) {
			return parameters.getParameterAt(selected).getValue();
		} else {
			return null;
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		String s = parameters.getParameterAt( 0 ).getName();//.getValue().toString();
		for( int i = 1; i < parameters.getNumberOfParameters(); i++ ) {
			s += ", " +  parameters.getParameterAt( i ).getName();//.getValue().toString();
		}
		return name + " (" + comment
			+ ", range={" + s + "}" 
			+ (defaultSelected>=0?", default = " + parameters.getParameterAt(defaultSelected).getValue():"")
			+ (isRequired() ? "" : ", OPTIONAL" )
			+ ")\t= " + getValue();
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		
		String selected = XMLParser.extractForTag( command, namePrefix ).toString();
		this.setValue( selected );
		if(this.getValue() instanceof GalaxyConvertible){
			((GalaxyConvertible)this.getValue()).fromGalaxy( namePrefix+"_opt"+getSelected(), command );
		}
	}
}