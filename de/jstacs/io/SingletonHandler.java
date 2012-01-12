package de.jstacs.io;

import java.lang.reflect.Field;

import de.jstacs.Singleton;

/**
 * 
 * @see de.jstacs.Singleton
 * 
 * @author Jens Keilwagen
 */
public class SingletonHandler {

	public static Singleton getSingelton( Class<? extends Singleton> singletonClass ) throws SecurityException, NoSuchFieldException, IllegalArgumentException, IllegalAccessException  {
		Field f = singletonClass.getField( "SINGLETON" );
		return (Singleton) (f.get(null));
	}
}
