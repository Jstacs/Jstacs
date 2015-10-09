package de.jstacs.utils;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.StringWriter;
import java.util.zip.DataFormatException;
import java.util.zip.DeflaterOutputStream;
import java.util.zip.InflaterInputStream;

import com.sun.org.apache.xml.internal.security.exceptions.Base64DecodingException;
import com.sun.org.apache.xml.internal.security.utils.Base64;


public class Compression {

	static{
		com.sun.org.apache.xml.internal.security.Init.init();
	}
	
	public static String unzip( String zipped ) throws DataFormatException, IOException, Base64DecodingException{
		
		InputStream in = new InflaterInputStream(new ByteArrayInputStream(Base64.decode(zipped)));
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
            byte[] buffer = new byte[8192];
            int len;
            while((len = in.read(buffer))>0){
                baos.write(buffer, 0, len);
            }
            return new String(baos.toByteArray(), "UTF-8");
		
	}
	
	public static String zip( String original ) throws IOException{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
        OutputStream out = new DeflaterOutputStream(baos);
        out.write(original.getBytes("UTF-8"));
        out.close();
        
		return Base64.encode(baos.toByteArray());
		
	}
	
	
	public static void main(String[] args) throws Exception {
		String text = "Hallo Welt! <tag> </tag>";
		
		String zipped = zip(text);
		
		System.out.println(zipped);
		
		String original = unzip( zipped );
		
		System.out.println(original);
		
	}
	
	
}
