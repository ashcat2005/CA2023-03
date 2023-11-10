
// hash numbers, prime numbers that ensure the hash is unique
int p1=73856093;
int p2=19349663;
int p3=83492791;

int m=10; // size of hash table. 
/*
where p1 = 73856093, p2 = 19349663, and p3 = 83492791 are
large prime numbers and where m is the hash table size. Please
note, that it generally cannot be avoided that several spatial cells
are mapped to the same hash value (hash collision). The effect of
overpopulated entries in the hash table might lead to a slow-down
of the neighborhood query. However, as suggested by Teschner the number of hash collisions can be reduced by
increasing the hash table size, i.e., trading memory for speed.
*/


