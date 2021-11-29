import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {
    public static ArrayList<String> readFasta() throws FileNotFoundException {

        ArrayList<String> database = new ArrayList<String>();
        try (Scanner sc = new Scanner(new File("src/tRNAs.fasta"))) {
            while (sc.hasNextLine()) {
                ArrayList<String> line = new ArrayList<String>();
                database.add(sc.nextLine().trim());
            }
        }
        return database;
    };
    public static void main(String[] args) throws FileNotFoundException {
        System.out.println(readFasta());
    }
}
