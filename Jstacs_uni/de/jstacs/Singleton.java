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
package de.jstacs;

import java.lang.reflect.Field;

/**
 * This interface states that the implementing class has only one immutable instance.
 * This instance has to be represented by the <code>public static final</code> field named <code>SINGLETON</code>.
 * No <code>public</code> constructor should be available.
 * 
 * {@link Singleton} instances can be saved and loaded using the {@link de.jstacs.io.XMLParser}.
 * 
 * @see de.jstacs.io.XMLParser
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface Singleton {
	
	/**
	 * This handler helps to retrieve the single instance of a {@link Singleton}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class SingletonHandler {
		
		/**
		 * This method helps to retrieve the single instance of a {@link Singleton} <code>singletonClass</code>.
		 * 
		 * @param singletonClass the class for which the 
		 * 
		 * @return the single instance of the given class
		 *  
		 * @throws SecurityException forwarded from {@link java.lang.Class#getField(String)}
		 * @throws NoSuchFieldException forwarded from {@link java.lang.Class#getField(String)}
		 * @throws IllegalArgumentException forwarded from {@link java.lang.reflect.Field#get(Object)} 
		 * @throws IllegalAccessException forwarded from {@link java.lang.reflect.Field#get(Object)}
		 */
		public static Singleton getSingelton( Class<? extends Singleton> singletonClass ) throws SecurityException, NoSuchFieldException, IllegalArgumentException, IllegalAccessException  {
			Field f = singletonClass.getField( "SINGLETON" );
			return (Singleton) (f.get(null));
		}
	}
	
}
