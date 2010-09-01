/* FIX THE BUG IN THE FOLLOWING FUNCTION, 
   CREATE A PATCH, AND SEND IT TO THE 
   MAILING LIST hermes2d@googlegroups.com.
   See the Sphinx tutorial for instructions. 
*/

int factorial(int n) {
  if(n == 0) return 0;
  else return n*factorial(n-1);
}
