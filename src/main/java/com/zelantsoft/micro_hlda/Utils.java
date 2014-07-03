package com.zelantsoft.micro_hlda;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Map;

import com.google.common.base.CharMatcher;

public class Utils {
	
	public static <T> void incInt(Map<T, Integer> map, T key)
	{
		Integer i = map.get(key);
		if(i == null)
		{
			i = new Integer(0);
		}
		i = new Integer(i.intValue() + 1);
		map.put(key, i);
	}
	
	public static <T> void addDouble(Map<T, Double> map, T key, double val)
	{
		Double i = map.get(key);
		if(i == null)
		{
			i = new Double(0);
		}
		i = new Double(i.doubleValue() + val);
		map.put(key, i);
	}
	
	public static void main(String [] args)
	{
		System.out.println(Utils.logb(532, 3.98));
	}
	
	public static String clean(String w)
	{
		return CharMatcher.JAVA_LETTER.retainFrom(w);
	}
	public static double logb( double a, double base )
	{
		return Math.log(a) / Math.log(base);
	}
	
	public static int teaseInt(String s)
	{
		StringBuilder sb = new StringBuilder();
		char [] chars = s.toCharArray();
		for(char ch : chars)
		{
			if(Character.isDigit(ch))
			{
				sb.append(ch);
			}	
		}
		if(sb.length() > 0)
		{
			return Integer.parseInt(sb.toString());
		}
		else
		{
			return 0;
		}
		
	}

	public static void downloadAsFile(InputStream is, File file) throws Exception
	{
		OutputStream fos = new FileOutputStream(file);
		download(is, fos);
		fos.close();
	}
	
	public static String downloadAsString(InputStream is) throws Exception
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		
		download(is, bos);
		
		return bos.toString("UTF-8");
	}
	
	public static void download(InputStream is, OutputStream os) throws Exception
	{
		int read = -1;
		byte [] buff = new byte[1024];
		
		while((read = is.read(buff)) > -1)
		{
			os.write(buff, 0, read);
		}
		os.flush();
	}
	
}
