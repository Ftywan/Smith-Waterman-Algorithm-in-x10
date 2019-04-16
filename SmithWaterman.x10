package smithwatermanalgo;
import x10.array.*;

public class SmityWaterman {
    private val seq1:String;
    private val seq2:String;
    private val length1:Int;
    private val length2:Int;
    private var score:Array_2[Double];
    private val MATCH_SCORE:Double;
    private val GAPOPEN_SCORE:Double;
    private val GAPEXT_SCORE:Double;
    private val DR_LEFT:Int = 1;// 0001
    private val DR_UP:Int = 2;// 0010
    private val DR_DIAG:Int = 4;// 0100
    private val DR_ZERO:Int = 8;// 1000
    private val BLOSUM62:Array_2[Long];
    private val MAP:HashMap[String, Int];
    private var sizesOfVerticalGaps:Array[Short];
    private var sizesOfHorizontalGaps:Array[Short];
    private var prevCells:Array_2[Double];
    

    public def SmithWatermans(seq1:Rail[String], seq2:Rail[String]) {
        this.seq1 = seq1;
        this.seq2 = seq2;
        this.length1 = seq1.length();
        this.length2 = seq2.length();
        score = new Array_2[Double](length1 + 1, length2 + 1);
        prevCells = new Array_2[Int](length1 + 1, length2 + 1);
   
    }

    public def similarity(i:Int, j:Int):Double {
        if (seq1.charAt(i - 1) == seq2.charAt(j - 1)) {
            return MATCH_SCORE;
        } else {
            return GAPOPEN_SCORE;
        }
    }

    public def initMatrix() {
        var i:Int = 0;
        var j:Int = 0;
        for (i in 0..(length1)) {
            score(i, 0) = 0;
            prevCells(i, 0) = DR_ZERO;
        }
        for (j in 0..(length2)) {
            score(0, j) = 0;
            prevCells(0, j) = DR_ZERO;
        }

    }
    public def readFastaFile(fastaFileName:String):String {
        val fastaFile = new File(fastaFileName);
        val fastaReader = new FileReader(fastaFile);
        val header = fastaReader.readLine();
        val builder = new StringBuilder();
        val line = fastaReader.readLine();
        while (line != null) {
            builder.append(line);
        } 
        return builder.toString();
    }

    public def readBlosumFile(blosumFileName:String):Array_2[Long] {
        val blosumFile = new File(blosumFileName);
        val blosumFileReader = blosumFile.openRead();
        val blosum62:Array_2[Long] = new Array_2[Long]()
        
      
    }

    public def initMap() {
        this.map = new HashMap[String, Int]();
        map.put("A", 0);
        map.put("R", 1);
        map.put("N", 2);
        map.put("D", 3);
        map.put("C", 4);
        map.put("Q", 5);
        map.put("E", 6);
        map.put("G", 7);
        map.put("H", 8);
        map.put("I", 9);
        map.put("L", 10);
        map.put("K", 11);
        map.put("M", 12);
        map.put("F", 13);
        map.put("P", 14);
        map.put("S", 15);
        map.put("T", 16);
        map.put("W", 17);
        map.put("X", 18);
        map.put("Y", 19);
        map.put("V", 20);
    }
}

/**
 * Blosum-62 substitution matrix
* #   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
* A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 
* R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 
* N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3 
* D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3 
* C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 
* Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2 
* E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2 
* G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 
* H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3 
* I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 
* L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 
* K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2 
* M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 
* F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 
* P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 
* S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2 
* T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 
* W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 
* Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 
* V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 
*/

