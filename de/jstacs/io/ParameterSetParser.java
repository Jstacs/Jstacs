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

package de.jstacs.io;

import java.lang.reflect.Constructor;

import de.jstacs.DataType;
import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;

/**
 * This class extracts values from {@link Parameter}s and creates instances of
 * {@link InstantiableFromParameterSet}s from a {@link ParameterSet}.
 * 
 * @author Jan Grau
 */
public class ParameterSetParser {

	/**
	 * Returns the <code>int</code> which is the value of the {@link Parameter}
	 * <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the <code>int</code> value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             <code>int</code>, i.e. its {@link DataType} is not
	 *             {@link DataType#INT}
	 * 
	 * @see DataType
	 * @see DataType#INT
	 */
	public static int getIntFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.INT ) {
			return (Integer)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type int." );
		}
	}

	/**
	 * Returns the <code>float</code> which is the value of the
	 * {@link Parameter} <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the <code>float</code> value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             <code>float</code>, i.e. its {@link DataType} is not
	 *             {@link DataType#FLOAT}
	 * 
	 * @see DataType
	 * @see DataType#FLOAT
	 */
	public static float getFloatFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.FLOAT ) {
			return (Float)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type float." );
		}
	}

	/**
	 * Returns the <code>double</code> which is the value of the
	 * {@link Parameter} <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the <code>double</code> value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             <code>double</code>, i.e. its {@link DataType} is not
	 *             {@link DataType#DOUBLE}
	 * 
	 * @see DataType
	 * @see DataType#DOUBLE
	 */
	public static double getDoubleFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.DOUBLE ) {
			return (Double)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type double." );
		}
	}

	/**
	 * Returns the <code>short</code> which is the value of the
	 * {@link Parameter} <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the <code>short</code> value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             <code>short</code>, i.e. its {@link DataType} is not
	 *             {@link DataType#SHORT}
	 * 
	 * @see DataType
	 * @see DataType#SHORT
	 */
	public static short getShortFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.SHORT ) {
			return (Short)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type short." );
		}
	}

	/**
	 * Returns the <code>long</code> which is the value of the {@link Parameter}
	 * <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the <code>long</code> value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             <code>long</code>, i.e. its {@link DataType} is not
	 *             {@link DataType#LONG}
	 * 
	 * @see DataType
	 * @see DataType#LONG
	 */
	public static long getLongFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.LONG ) {
			return (Long)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type long." );
		}
	}

	/**
	 * Returns the <code>byte</code> which is the value of the {@link Parameter}
	 * <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the <code>byte</code> value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             <code>byte</code>, i.e. its {@link DataType} is not
	 *             {@link DataType#BYTE}
	 * 
	 * @see DataType
	 * @see DataType#BYTE
	 */
	public static byte getByteFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.BYTE ) {
			return (Byte)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type byte." );
		}
	}

	/**
	 * Returns the <code>boolean</code> which is the value of the
	 * {@link Parameter} <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the <code>boolean</code> value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             <code>boolean</code>, i.e. its {@link DataType} is not
	 *             {@link DataType#BOOLEAN}
	 * 
	 * @see DataType
	 * @see DataType#BOOLEAN
	 */
	public static boolean getBooleanFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.BOOLEAN ) {
			return (Boolean)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type boolean." );
		}
	}

	/**
	 * Returns the {@link String} which is the value of the {@link Parameter}
	 * <code>par</code>.
	 * 
	 * @param par
	 *            the {@link Parameter}
	 * 
	 * @return the {@link String} value of the {@link Parameter}
	 * 
	 * @throws WrongParameterTypeException
	 *             if <code>par</code> is not a {@link Parameter} of type
	 *             {@link String}, i.e. its {@link DataType} is not
	 *             {@link DataType#STRING}
	 * 
	 * @see DataType
	 * @see DataType#STRING
	 */
	public static String getStringFromParameter( Parameter par ) throws WrongParameterTypeException {
		if( par.getDatatype() == DataType.STRING ) {
			return (String)par.getValue();
		} else {
			throw new WrongParameterTypeException( "Parameter " + par.getName() + " is not of type String." );
		}
	}

	/**
	 * Returns an instance of a subclass of {@link InstantiableFromParameterSet}
	 * that can be instantiated by the {@link InstanceParameterSet}
	 * <code>pars</code>. The instance class is taken from <code>pars</code> via
	 * the method {@link InstanceParameterSet#getInstanceClass()}.
	 * 
	 * @param pars
	 *            the {@link InstanceParameterSet}
	 * 
	 * @return the instance
	 * 
	 * @throws NotInstantiableException
	 *             if {@link InstanceParameterSet#getInstanceClass()} of
	 *             <code>pars</code> is <code>null</code>, could not be found or
	 *             cannot be instantiated from <code>pars</code>
	 * 
	 * @see InstanceParameterSet#getInstanceClass()
	 * @see ParameterSetParser#getInstanceFromParameterSet(ParameterSet, Class)
	 */
	public static InstantiableFromParameterSet getInstanceFromParameterSet( InstanceParameterSet pars ) throws NotInstantiableException {
		return getInstanceFromParameterSet( pars, pars.getInstanceClass() );
	}

	/**
	 * Returns an instance of a subclass of {@link InstantiableFromParameterSet}
	 * that can be instantiated by the {@link ParameterSet} <code>pars</code>.
	 * The instance class is taken from <code>instanceClass</code>.
	 * 
	 * @param pars
	 *            the {@link ParameterSet}
	 * @param instanceClass
	 *            the class that shall be instantiated
	 * 
	 * @return the instance
	 * 
	 * @throws NotInstantiableException
	 *             if <code>instanceClass</code> is <code>null</code>, could not
	 *             be found or cannot be instantiated from <code>pars</code>
	 */
	public static InstantiableFromParameterSet getInstanceFromParameterSet( ParameterSet pars, Class instanceClass ) throws NotInstantiableException {
		if( instanceClass == null ) {
			throw new NotInstantiableException( "An instance class must be specified." );
		}
		Class parClass = pars.getClass();
		Constructor construct = null;
		while( construct == null && parClass != Object.class ) {
			try {
				construct = instanceClass.getConstructor( new Class[]{ parClass } );
			} catch ( Exception e ) {
				construct = null;
				parClass = parClass.getSuperclass();
			}
		}
		if( construct == null ) {
			throw new NotInstantiableException( "No appropriate constructor found." );
		} else {
			try {
				return (InstantiableFromParameterSet)construct.newInstance( new Object[]{ pars } );
			} catch ( Exception e ) {
				throw new NotInstantiableException( e.getCause().getMessage() );
			}
		}
	}

	/**
	 * An {@link Exception} that is thrown if an instance of some class could
	 * not be created.
	 * 
	 * @author Jan Grau
	 */
	public static class NotInstantiableException extends ParameterException {

		private static final long serialVersionUID = 1L;

		/**
		 * Creates a new instance of a {@link NotInstantiableException} with a
		 * given error message.
		 * 
		 * @param msg
		 *            the error message
		 * 
		 * @see ParameterException#ParameterException(String)
		 */
		public NotInstantiableException( String msg ) {
			super( msg );
		}

	}

	/**
	 * An {@link Exception} that is thrown if the {@link DataType} of a
	 * {@link Parameter} is not appropriate for some purpose.
	 * 
	 * @author Jan Grau
	 * 
	 * @see de.jstacs.parameters.Parameter
	 * @see DataType
	 */
	public static class WrongParameterTypeException extends ParameterException {

		private static final long serialVersionUID = 1L;

		/**
		 * Creates a new instance of a {@link WrongParameterTypeException} with
		 * a given error message.
		 * 
		 * @param msg
		 *            the error message
		 * 
		 * @see ParameterException#ParameterException(String)
		 */
		public WrongParameterTypeException( String msg ) {
			super( msg );
		}
	}
}
