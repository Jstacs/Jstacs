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

package de.jstacs.utils;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.zip.DeflaterOutputStream;
import java.util.zip.InflaterInputStream;

//import javax.xml.bind.DatatypeConverter;
import java.util.Base64;

/**
 * Class for compressing and de-compressing {@link String}s
 * using ZIP.
 * @author Jan Grau
 *
 */
public class Compression {

	/**
	 * De-compressed the original {@link String} from the supplied (Base64) compressed {@link String}.
	 * @param zipped the zipped {@link String}
	 * @return the de-compressed {@link String}
	 * @throws IOException if the {@link String} could not be de-compressed
	 */
	public static String unzip( String zipped ) throws IOException {
		
		
		//InputStream in = new InflaterInputStream(new ByteArrayInputStream(Base64.decode(zipped)));
		//InputStream in = new InflaterInputStream(new ByteArrayInputStream(DatatypeConverter.parseBase64Binary(zipped)));
		InputStream in = new InflaterInputStream(new ByteArrayInputStream(Base64.getDecoder().decode(zipped)));
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
            byte[] buffer = new byte[8192];
            int len;
            while((len = in.read(buffer))>0){
                baos.write(buffer, 0, len);
            }
            return new String(baos.toByteArray(), "UTF-8");
		
	}
	
	/**
	 * Compresses the supplied original {@link String}.
	 * @param original the original {@link String}
	 * @return the compressed and Base64 encoded {@link String}
	 * @throws IOException if the {@link String} could not be compressed
	 */
	public static String zip( String original ) throws IOException{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
        OutputStream out = new DeflaterOutputStream(baos);
        out.write(original.getBytes("UTF-8"));
        out.close();
        
		//return Base64.encode(baos.toByteArray());
        //return DatatypeConverter.printBase64Binary(baos.toByteArray());
		return Base64.getEncoder().encodeToString(baos.toByteArray());
	}
	
	
}
