package de.jstacs.tools;

import java.io.IOException;
import java.io.OutputStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ProtocolOutputStream extends OutputStream {

	private Protocol prot;
	boolean toWarning;

	StringBuffer buffer;
	Pattern pat;
	
	public ProtocolOutputStream(Protocol prot, boolean toWarning){
		this(prot,toWarning,"([\\n\\s]|.)*");
	}
	
	public ProtocolOutputStream(Protocol prot, boolean toWarning, String regExp) {
		this.prot = prot;
		this.toWarning = toWarning;
		buffer = new StringBuffer();
		pat = Pattern.compile(regExp);
	}
	
	@Override
	public void write(int b) throws IOException {
		buffer.append(new String(new byte[]{(byte)b}));
		
		Matcher m = pat.matcher(buffer);
		int end = 0;
		while(m.find()){
			String gr = m.group(0); 
			if(toWarning){
				prot.appendWarning(gr);
			}else{
				prot.append(gr);
			}
			end = m.end();
		}
		buffer.delete(0, end);
	}

	@Override
	public void write(byte[] b, int off, int len) throws IOException {
		buffer.append(new String(b,off,len));

		Matcher m = pat.matcher(buffer);
		int end = 0;
		while(m.find()){
			String gr = m.group(0); 
			if(toWarning){
				prot.appendWarning(gr);
			}else{
				prot.append(gr);
			}
			end = m.end();
		}
		buffer.delete(0, end);
	}

	@Override
	public void flush() throws IOException {

	}

	@Override
	public void close() throws IOException {

	}
	
	
	

}
