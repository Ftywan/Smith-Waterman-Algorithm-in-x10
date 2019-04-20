import x10.lang.String;
import x10.io.Console;
import x10.io.File;
import x10.util.StringBuilder;
import x10.io.FileReader;
import x10.lang.Exception;


public class FastaReader {

    public static def main(argv: Rail[String]) {
        val reader = new FastaReader();

        val content:String = reader.readFastaFile(argv(0));
        Console.OUT.println(content);
    }
    public def readFastaFile(fastaFileName:String):String {
        val fastaFile:File = new File(fastaFileName);
        val fastaReader:FileReader = fastaFile.openRead();
        val header:String = fastaReader.readLine();
        val builder:StringBuilder = new StringBuilder();
        var line:String = null;
        try {
            line = fastaReader.readLine().trim();
           
        } catch (e:Exception) {
            Console.OUT.println("Exception has been caught");
        }

        while (true) {
            builder.add(line);
            try {
                line = fastaReader.readLine().trim();
            } catch (e:Exception) {
                Console.OUT.println("Exception has been caught");
                break;
            }
        }
        return builder.toString();
    }
}
