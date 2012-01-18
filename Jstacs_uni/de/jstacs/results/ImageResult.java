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

package de.jstacs.results;

import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import javax.imageio.ImageIO;

import com.sun.org.apache.xml.internal.security.utils.Base64;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * A class for results that are images of the PNG format. The images themselves
 * cannot be stored to an XML representation and thus only the description
 * (name, etc.) is stored.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ImageResult extends Result {

	/**
	 * The image as {@link BufferedImage}
	 */
	private BufferedImage image;

	/**
	 * Constructs a new {@link ImageResult} from a {@link BufferedImage}.
	 * 
	 * @param name
	 *            the name of the image
	 * @param comment
	 *            a comment on the image
	 * @param image
	 *            the image itself
	 */
	public ImageResult(String name, String comment, BufferedImage image) {
		super( name, comment, DataType.PNG );
		this.image = image;
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ImageResult} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ImageResult} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ImageResult(StringBuffer xml) throws NonParsableException {
		super( xml );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.Result#getResult()
	 */
	@Override
	public BufferedImage getValue() {
		return image;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "ImageResult";
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer xml ) {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(1000);
		try {
			ImageIO.write( image, "png", baos );
		} catch ( IOException e ) {
			throw new RuntimeException( e.getMessage() );
		}
		XMLParser.appendObjectWithTags(xml, Base64.encode(baos.toByteArray()), "image");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		try {
			byte[] bytearray = Base64.decode( (String) XMLParser.extractObjectForTags( representation, "image" ) );
			image = ImageIO.read( new ByteArrayInputStream( bytearray ) );
		} catch (Exception e) {
			throw new NonParsableException( e.getMessage() );
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return name + ": [image]";
	}
}