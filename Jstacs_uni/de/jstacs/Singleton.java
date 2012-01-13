package de.jstacs;

/**
 * This interface states that the implementing class has only one immutable instance.
 * This instance has to be represented by the <code>public static final</code> field named <code>SINGLETON</code>.
 * No <code>public</code> constructor should be available.
 * 
 * {@link Singleton} instances can be saved and loaded using the {@link de.jstacs.io.XMLParser}.
 * 
 * @see de.jstacs.io.XMLParser
 * @see de.jstacs.io.SingletonHandler
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface Singleton {

}
