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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.io.XMLParser;

/**
 * (Container) class for a set of {@link Parameter}s. This is the base class for
 * other {@link ParameterSet}s that provide more specialized methods, e.g. for
 * expanding the set of {@link Parameter}s at runtime (
 * {@link ExpandableParameterSet} ) or defining an array of {@link ParameterSet}
 * s from a common template ( {@link ArrayParameterSet} ).
 * 
 * @author Jan Grau
 */
public abstract class ParameterSet implements Storable, Cloneable,
		RangeIterator, GalaxyConvertible {
	
	/**
	 * This method tries to find the correct name ({@link String}) for your
	 * choice. This method is useful if you handle {@link CollectionParameter}s.
	 * 
	 * @param names
	 *            the names
	 * @param values
	 *            the values that can be set
	 * @param current
	 *            the value to be set
	 * @param hasAlternative
	 *            indicates whether the last entry of names is an alternative
	 *            parameter
	 * 
	 * @return the {@link StringIndexOutOfBoundsException} that has to be used
	 * 
	 * @throws IllegalArgumentException
	 *             if no match could be found
	 * 
	 * @see CollectionParameter
	 */
	public static int getIndex(String[] names, Object[] values,
			Comparable current, boolean hasAlternative)
			throws IllegalArgumentException {
		int i = 0, l = values.length;
		if (hasAlternative) {
			l--;
		}
		while (i < l && !current.equals(values[i])) {
			i++;
		}
		if (i == values.length) {
			throw new IllegalArgumentException(
					"Could not find a matching constant.");
		} else {
			return i;
		}
	}
	
	/**
	 * The set of parameters
	 */
	protected ParameterList parameters;

	/**
	 * The error message of the last error or <code>null</code>
	 */
	protected String errorMessage;

	/**
	 * Indicates if the {@link Parameter}s of this {@link ParameterSet} that
	 * implement {@link Rangeable} and return <code>true</code> shall be
	 * replaced by their ranged instances
	 */
	protected boolean ranged = false;

	/**
	 * If this {@link ParameterSet} is contained in a
	 * {@link ParameterSetContainer}, this variable holds a reference to that
	 * {@link ParameterSetContainer}.
	 */
	protected ParameterSetContainer parent;

	private long id;

	/**
	 * Constructs a new {@link ParameterSet} with empty parameter values. The
	 * set of parameters is loaded by the method {@link #loadParameters()} as
	 * soon as one of the accession methods, e.g.
	 * {@link #getNumberOfParameters()}, is called.
	 */
	public ParameterSet() {
		this.id = System.currentTimeMillis() + this.hashCode();
	}

	/**
	 * Creates a full clone (deep copy) of this {@link ParameterSet}. As a
	 * convenience-method the user can use
	 * <code>fillWithStandardFieldsForClone(ParameterSet)</code> on a newly
	 * created instance of a subclass of {@link ParameterSet} to obtain a
	 * clone/copy of all standard member variables (those already defined in
	 * {@link ParameterSet}) in the passed {@link ParameterSet}. Using this
	 * method, the cloning process becomes merely three-step:<br>
	 * <ul>
	 * <li>Create a new instance of your subclass of {@link ParameterSet}, most
	 * likely with an empty constructor or the one taking just the instance
	 * class.
	 * <li>Call <code>this.fillWithStandardFieldsForClone</code> on this
	 * instance.
	 * <li>Return the instance.
	 * </ul>
	 * This method fulfills the conventions of {@link Object}'s method
	 * {@link Object#clone()}.
	 * 
	 * @return a deep clone/copy of this {@link ParameterSet}
	 */
	@Override
	public ParameterSet clone() throws CloneNotSupportedException {

		ParameterSet ret = (ParameterSet) super.clone();

		if (this.parameters != null) {
			ret.initParameterList(this.parameters.size());
			for (int i = 0; i < this.parameters.size(); i++) {
				ret.parameters.add(i, this.parameters.get(i).clone());
			}
			Iterator<Parameter> it = ret.parameters.iterator();
			while (it.hasNext()) {
				Parameter par = it.next();
				if (par.getNeededReferenceId() != null) {
					ParameterSet set = findParameterSet(par
							.getNeededReferenceId(), false);
					if (set != null) {
						par.setNeededReference(set);
					}
				}
			}
		}
		propagateId(this.id, this, false);
		return ret;
	}

	/**
	 * Initializes the internal set of {@link Parameter}s, which is a
	 * {@link ParameterList}.
	 */
	protected final void initParameterList() {
		this.parameters = new ParameterList();
	}

	/**
	 * Initializes the internal set of {@link Parameter}s, which is a
	 * {@link ParameterList}, with an initial number of {@link Parameter}s of
	 * <code>initCapacity</code>.
	 * 
	 * @param initCapacity
	 *            the initial number of {@link Parameter}s
	 */
	protected final void initParameterList(int initCapacity) {
		this.parameters = new ParameterList(initCapacity);
	}

	/**
	 * Loads the parameters for this {@link ParameterSet}. This is in most cases
	 * done by calling {@link #initParameterList()} or
	 * {@link #initParameterList(int)} to initialize
	 * {@link ParameterSet#parameters} and afterwards filling
	 * {@link ParameterSet#parameters} with instances of subclasses of
	 * {@link Parameter}.
	 * 
	 * @throws Exception
	 *             if the parameters could not be loaded
	 * 
	 * @see Parameter
	 * @see ParameterSet#parameters
	 * @see ParameterSet#initParameterList()
	 * @see ParameterSet#initParameterList(int)
	 */
	protected abstract void loadParameters() throws Exception;

	/**
	 * Returns <code>true</code> if the parameters of this {@link ParameterSet}
	 * have already been loaded using {@link #loadParameters()}.
	 * 
	 * @return the state of the parameters (loaded or not)
	 */
	public boolean parametersLoaded() {
		return (parameters != null);
	}

	/**
	 * Returns the index of the first parameter in the set of required
	 * parameters that has not been set. If all required parameters have been
	 * set, it returns -1.
	 * 
	 * @return the index of the required parameter that has not been set
	 */
	/*
	 * public int getIndexOfRequiredParameterNotSet() { if( requiredParameters
	 * == null ) { return -1; } else { for( int i = 0; i <
	 * requiredParameters.length; i++ ) { if( !requiredParameters[i].isSet() ) {
	 * return i; } } return -1; } }
	 */

	/**
	 * Returns <code>true</code> if this {@link ParameterSet} contains only
	 * atomic parameters, i.e. the parameters do not contain
	 * {@link ParameterSet}s themselves.
	 * 
	 * @return <code>true</code> if the {@link ParameterSet} contains only
	 *         atomic parameters, <code>false</code> otherwise
	 */
	public boolean isAtomic() {
		if (parameters == null) {
			try {
				loadParameters();
				if (ranged) {
					replaceParametersWithRangedInstance();
				}
			} catch (Exception e) {
				e.printStackTrace();
				return false;
			}
		}
		for (int i = 0; i < parameters.size(); i++) {
			if (!parameters.get(i).isAtomic()) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns <code>true</code> if all parameters in this {@link ParameterSet}
	 * are either set by the user or have default values. If any additional
	 * constraints are required on your parameters you should either specify
	 * some {@link de.jstacs.parameters.validation.ParameterValidator} on these parameters or implement these
	 * constraints by overriding this method in your implementation of
	 * {@link ParameterSet}. It is recommended to specify a useful remark which
	 * constraint failed in the member-variable <code>errorMessage</code>, which
	 * will be displayed to the user. In the overriding method
	 * <code>super.</code>{@link #hasDefaultOrIsSet()} should be called prior to
	 * checking specific constraints.
	 * 
	 * @return <code>true</code> if all parameters have some allowed value set,
	 *         <code>false</code> otherwise
	 */
	public boolean hasDefaultOrIsSet() {
		if (parameters == null) {
			try {
				loadParameters();
				if (ranged) {
					replaceParametersWithRangedInstance();
				}
			} catch (Exception e) {
				errorMessage = "Parameters could not be loaded";
				e.printStackTrace();
				return false;
			}
		}
		for (int i = 0; i < parameters.size(); i++) {
			if (parameters.get(i).isRequired()
					&& (!parameters.get(i).hasDefaultOrIsSet())) {
				// errorMessage =
				// "Parameter no."+i+" ("+parameters.get(i).getName()+") was not set.";
				/*
				 * if(parameters.get(i).getErrorMessage() != null){ errorMessage
				 * =
				 * "Parameter no. "+(i+1)+" had the following error: "+parameters
				 * .get(i).getErrorMessage(); }
				 */
				return false;
			}
		}
		errorMessage = null;
		return true;
	}

	/**
	 * Returns the message of the last error that occurred. If no error occurred
	 * this method returns <code>null</code>.
	 * 
	 * @return the error message of the last error or <code>null</code>
	 */
	public String getErrorMessage() {
		return errorMessage;
	}

	/**
	 * Constructs a {@link ParameterSet} out of an array of {@link Parameter}s.
	 * The {@link Parameter}s are not cloned, but passed by reference.
	 * 
	 * @param parameters
	 *            the {@link Parameter}s
	 */
	protected ParameterSet(Parameter[] parameters) {
		this.id = System.currentTimeMillis() + this.hashCode();
		initParameterList(parameters.length);
		for (int i = 0; i < parameters.length; i++) {
			this.parameters.add(parameters[i]);
		}
	}

	/**
	 * Constructs a {@link ParameterSet} out of an {@link ArrayList} of
	 * {@link Parameter}s. The {@link Parameter}s are not cloned, but passed by
	 * reference.
	 * 
	 * @param parameters
	 *            the {@link Parameter}s
	 */
	protected ParameterSet(ArrayList<Parameter> parameters) {
		this.id = System.currentTimeMillis() + this.hashCode();
		initParameterList();
		Iterator<Parameter> parIt = parameters.iterator();
		while (parIt.hasNext()) {
			this.parameters.add(parIt.next());
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ParameterSet} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>representation</code>
	 */
	public ParameterSet(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/**
	 * Returns the number of parameters in the {@link ParameterSet}.
	 * 
	 * @return the number of parameters in the {@link ParameterSet}
	 */
	public int getNumberOfParameters() {
		if (parameters == null) {
			try {
				loadParameters();
				if (ranged) {
					replaceParametersWithRangedInstance();
				}
			} catch (Exception e) {
				e.printStackTrace();
				return 0;
			}
		}
		if (parameters != null) {
			return parameters.size();
		} else {
			return 0;
		}
	}

	/**
	 * Returns the {@link Parameter} at position <code>i</code>.
	 * 
	 * @param i
	 *            the position in the {@link ParameterList}
	 * 
	 * @return the {@link Parameter} at position <code>i</code>
	 */
	public Parameter getParameterAt(int i) {
		if (parameters == null) {
			try {
				loadParameters();
				if (ranged) {
					replaceParametersWithRangedInstance();
				}
			} catch (Exception e) {
				e.printStackTrace();
				return null;
			}
		}
		return parameters.get(i);
	}

	/**
	 * Replaces all {@link Parameter}s in this {@link ParameterSet} by their
	 * equivalents implementing the {@link Rangeable} interface.
	 * 
	 * @throws Exception
	 *             if these instances could not be created
	 * 
	 * @see ParameterSet#replaceParametersWithRangedInstance()
	 */
	public void makeRanged() throws Exception {
		if (parameters != null) {
			replaceParametersWithRangedInstance();
		}
		this.ranged = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#next()
	 */
	public boolean next() throws ParameterException {
		if (ranged) {
			RangeIterator ri;
			for (int i = 0; i < parameters.size(); i++) {
				if (parameters.get(i) instanceof RangeIterator) {
					ri = (RangeIterator) parameters.get(i);
					if (ri.next()) {
						return true;
					} else {
						ri.resetToFirst();
					}
				}
			}
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#resetToFirst()
	 */
	public void resetToFirst() {
		if (ranged) {
			for (int i = 0; i < parameters.size(); i++) {
				if (parameters.get(i) instanceof RangeIterator) {
					((RangeIterator) parameters.get(i)).resetToFirst();
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#getNumberOfValues()
	 */
	public int getNumberOfValues() {
		if (ranged) {
			int erg = 1;
			for (int i = 0; i < parameters.size(); i++) {
				if (parameters.get(i) instanceof RangeIterator) {
					erg *= ((RangeIterator) parameters.get(i))
							.getNumberOfValues();
				}
			}
			return erg;
		} else {
			return 1;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#isRanged()
	 */
	public boolean isRanged() {
		if (ranged) {
			return getNumberOfValues() > 1;
		} else {
			return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#valuesToString()
	 */
	public String valuesToString() {
		String erg = null;
		for (int i = 0; i < parameters.size(); i++) {
			if (ranged && parameters.get(i) instanceof RangeIterator) {
				if (erg == null) {
					erg = ((RangeIterator) parameters.get(i)).valuesToString();
				} else {
					erg += " x "
							+ ((RangeIterator) parameters.get(i))
									.valuesToString();
				}
			} else {
				if (erg == null) {
					erg = "[" + parameters.get(i).getValue().toString() + "]";
				} else {
					erg += " x [" + parameters.get(i).getValue().toString()
							+ "]";
				}
			}
		}
		return erg;
	}

	/**
	 * Replaces all {@link Parameter}s in this {@link ParameterSet} by their
	 * equivalents implementing the {@link Rangeable} interface.
	 * 
	 * @throws Exception
	 *             if these instances could not be created
	 * 
	 * @see ParameterSet#makeRanged()
	 */
	protected void replaceParametersWithRangedInstance() throws Exception {
		for (int i = 0; i < parameters.size(); i++) {
			if (parameters.get(i) instanceof Rangeable) {
				Rangeable r = (Rangeable) parameters.get(i);
				if (r.isRangeable()) {
					parameters.set(i, r.getRangedInstance());
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		if (parameters != null) {
			StringBuffer buf2 = new StringBuffer();
			XMLParser.appendObjectWithTags(buf2, parameters.size(), "numberOfParameters");
			for (int i = 0; i < parameters.size(); i++) {
				XMLParser.appendObjectWithTags(buf2, parameters.get( i ), "parameter");
				/*
				if (parameters.get(i) == null) {
					XMLParser.appendObjectWithTags(buf2, "null", "parameter");
				} else {
					StringBuffer buf3 = new StringBuffer();
					XMLParser.appendObjectWithTags(buf3, parameters.get(i).getClass().getName(), "className");
					buf3.append(parameters.get(i).toXML());
					XMLParser.addTags(buf3, "parameter");
					buf2.append(buf3);
				}
				*/
			}
			XMLParser.addTags(buf2, "set");
			buf.append(buf2);
		} else {
			XMLParser.appendObjectWithTags(buf, null, "set");
		}
		XMLParser.appendObjectWithTags(buf, ranged, "ranged");
		XMLParser.appendObjectWithTags(buf, id, "id");
		XMLParser.addTags(buf, "parameterSet");

		return buf;
	}

	/**
	 * Simplifies all {@link Parameter}s in this {@link ParameterSet}.
	 * 
	 * @see Parameter#simplify()
	 */
	public void simplify() {
		if (parameters == null) {
			return;
		}
		for (int i = 0; i < parameters.size(); i++) {
			if (parameters.get(i) instanceof CollectionParameter) {
				parameters.get(i).simplify();
			} else if (parameters.get(i).getValue() instanceof ParameterSet) {
				((ParameterSet) parameters.get(i).getValue()).simplify();
			}
		}
	}

	/**
	 * Resets all {@link Parameter}s in this {@link ParameterSet} to their
	 * default values or <code>null</code> if not default value was provided.
	 * 
	 * @see Parameter#reset()
	 */
	public void reset() {
		/*
		 * if( parameters == null ) { return; } for( int i = 0; i <
		 * parameters.size(); i++ ) { parameters.get( i ).reset(); }
		 */
		parameters = null;
	}

	/**
	 * Parses the instance fields of a {@link ParameterSet} from the XML
	 * representation as returned by {@link #toXML()}.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 * 
	 * @see ParameterSet#toXML()
	 */
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		
		representation = XMLParser.extractForTag(representation, "parameterSet");
		StringBuffer buf = XMLParser.extractForTag(representation, "set");
		if (!buf.toString().equalsIgnoreCase("null")) {
			int numPars = XMLParser.extractObjectForTags(buf, "numberOfParameters", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			parameters = new ParameterList(numPars);
			for (int i = 0; i < numPars; i++) {
				parameters.add( (Parameter) XMLParser.extractObjectForTags( buf, "parameter" ) );

				/*
				String s = buf.toString();
				String par = XMLParser.extractObjectAndAttributesForTags(buf, "parameter", null,null,String.class,false );
				if (par.equals("null")) {
					parameters.add(null);
				} else {
					if( buf2 == null ) {
						buf2 = new StringBuffer(par);
					} else {
						buf2.delete(0, buf2.length());
						buf2.append(par);
					}
					String className = XMLParser.extractObjectForTags(buf2, "className", String.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
					s = buf2.toString();
					try {
						parameters.add((Parameter) Class.forName(className)
								.getConstructor( new Class[] { StringBuffer.class }).newInstance(buf2));
					} catch (Exception e) {
						System.out.println();
						System.out.println( "StringBuffer:\n" + s );
						System.out.println();
						
						System.out.println("parameter: " + i);
						System.out.println("class    : " + className);
						System.out.println("buffer   : " + par);
						
						System.out.println();
						System.out.println("buffer remaining:\n" + buf2 );
						
						try {
							Thread.sleep(50);
						} catch (InterruptedException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
						
						NonParsableException n = new NonParsableException(e.getMessage());
						n.setStackTrace(e.getStackTrace());
						System.out.print("ParameterSet: ");
						e.printStackTrace();
						
						System.exit(1);
						throw n;
					}
				}
				/**/
			}
		}
		StringBuffer help = XMLParser.extractForTag(representation, "ranged");
		if (help == null) {
			ranged = false;
		} else {
			ranged = Boolean.parseBoolean(help.toString());
		}
		id = XMLParser.extractObjectForTags(representation, "id", long.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		propagateId(id, this, false);
		recieveId();
	}

	/**
	 * Searches for all {@link Parameter#neededReferenceId}s in the hierarchy
	 * below this {@link ParameterSet} and sets the corresponding
	 * {@link Parameter#neededReference}s.
	 */
	protected void recieveId() {
		if (parameters != null) {
			Iterator<Parameter> it = parameters.iterator();
			while (it.hasNext()) {
				Parameter par = it.next();
				Long neededId = par.getNeededReferenceId();
				if (neededId != null) {
					par.setNeededReference(findParameterSet(neededId, false));
				}
			}
		}
	}

	private ParameterSet findParameterSet(long id, boolean fwd) {
		if (id == this.id) {
			return this;
		} else {
			if (!fwd && this.getParent() != null
					&& this.getParent().getParent() != null) {
				ParameterSet set = this.getParent().getParent()
						.findParameterSet(id, false);
				if (set != null) {
					return set;
				} else {
					return null;
				}
			} else if (parameters != null) {
				Iterator<Parameter> it = parameters.iterator();
				while (it.hasNext()) {
					Parameter par = it.next();
					if (par instanceof ParameterSetContainer) {
						ParameterSet set = ((ParameterSetContainer) par)
								.getValue().findParameterSet(id, true);
						if (set != null) {
							return set;
						}
					}
				}
			}

			return null;
		}
	}

	/**
	 * Propagates the id of this {@link ParameterSet} to all {@link Parameter}s
	 * above and below in the hierarchy.
	 */
	protected void propagateId() {
		propagateId(id, this, false);
	}

	private boolean propagateId(long id, ParameterSet set, boolean fwd) {
		if (!fwd && this.getParent() != null
				&& this.getParent().getParent() != null) {
			return this.getParent().getParent().propagateId(id, set, false);
		} else if (parameters != null) {
			Iterator<Parameter> it = this.parameters.iterator();
			while (it.hasNext()) {
				Parameter par = it.next();
				if (par.getNeededReferenceId() != null
						&& par.getNeededReferenceId() == id) {
					par.setNeededReference(set);
					return true;
				}
				if (par instanceof ParameterSetContainer) {
					boolean re = ((ParameterSetContainer) par).getValue()
							.propagateId(id, set, true);
					if (re) {
						return true;
					}
				}
			}
			return false;
		}
		return false;
	}

	/**
	 * Returns the id of this {@link ParameterSet}.
	 * 
	 * @return the id of the {@link ParameterSet}
	 */
	public long getId() {
		return id;
	}

	/**
	 * Returns the enclosing {@link ParameterSetContainer} of this
	 * {@link ParameterSet} or <code>null</code> if none exists.
	 * 
	 * @return the enclosing {@link ParameterSetContainer} (<code>parent</code>)
	 */
	public ParameterSetContainer getParent() {
		return parent;
	}

	/**
	 * Sets the enclosing {@link ParameterSetContainer} of this
	 * {@link ParameterSet} to <code>parent</code>.
	 * 
	 * @param parent
	 *            the new parent
	 */
	public void setParent(ParameterSetContainer parent) {
		this.parent = parent;
	}

	/**
	 * Class for a {@link java.util.List} of {@link Parameter}s that basically
	 * has the same functionality as {@link ArrayList}, but additionally takes
	 * care of the references {@link Parameter#parent}.
	 * 
	 * @author Jan Grau
	 */
	public class ParameterList extends ArrayList<Parameter> {

		private static final long serialVersionUID = 5646161850340681495L;

		/**
		 * Creates a new empty {@link ParameterList}.
		 */
		public ParameterList() {
			super();
		}

		/**
		 * Creates a new {@link ParameterList} from an existing
		 * {@link Collection} of {@link Parameter}s. This may be another
		 * {@link ParameterList}.
		 * 
		 * @param c
		 *            the {@link Collection} of {@link Parameter}s
		 */
		public ParameterList(Collection<? extends Parameter> c) {
			super(c);
			Iterator<? extends Parameter> it = c.iterator();
			while (it.hasNext()) {
				it.next().setParent(ParameterSet.this);
			}
		}

		/**
		 * Creates a new empty {@link ParameterList} with a defined initial
		 * capacity.
		 * 
		 * @param initialCapacity
		 *            the initial capacity of the new {@link ParameterList}
		 */
		public ParameterList(int initialCapacity) {
			super(initialCapacity);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.ArrayList#add(int, java.lang.Object)
		 */
		@Override
		public void add(int index, Parameter element) {
			element.setParent(ParameterSet.this);
			super.add(index, element);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.ArrayList#add(java.lang.Object)
		 */
		@Override
		public boolean add(Parameter o) {
			o.setParent(ParameterSet.this);
			return super.add(o);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.ArrayList#addAll(java.util.Collection)
		 */
		@Override
		public boolean addAll(Collection<? extends Parameter> c) {
			Iterator<? extends Parameter> it = c.iterator();
			while (it.hasNext()) {
				it.next().setParent(ParameterSet.this);
			}
			return super.addAll(c);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.ArrayList#addAll(int, java.util.Collection)
		 */
		@Override
		public boolean addAll(int index, Collection<? extends Parameter> c) {
			Iterator<? extends Parameter> it = c.iterator();
			while (it.hasNext()) {
				it.next().setParent(ParameterSet.this);
			}
			return super.addAll(index, c);
		}

	}
	
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer ) throws Exception {
		StringBuffer xml = new StringBuffer();
		for(int i=0;i<getNumberOfParameters();i++){
			((GalaxyConvertible)getParameterAt( i )).toGalaxy( namePrefix+"_ps", configPrefix, depth+1, descBuffer, configBuffer );
			descBuffer.append( "\n" );
			configBuffer.append( "\n" );
		}	
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		for(int i=0;i<getNumberOfParameters();i++){
			((GalaxyConvertible)getParameterAt( i )).fromGalaxy( namePrefix+"_ps", command );
		}
	}
	
	

}
