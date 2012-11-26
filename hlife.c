/*
 * hlife 0.96 by Radical Eye Software.
 *
 * Usage:
 *   hlife [-i increment] [-m maxgen] [-r rules] [-o outfile]
 *         [-M maxmem] [-q] [-2] lifefile
 *
 *   All good ideas here were originated by Gosper or Bell or others, I'm
 *   sure, and all bad ones by yours truly.
 *
 *   The main reason I wrote this program was to attempt to push out the
 *   evaluation of metacatacryst as far as I could.  So this program
 *   really does very little other than compute life as far into the
 *   future as possible, using as little memory as possible (and reusing
 *   it if necessary).  No UI, few options.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
/*
 *   Into instances of this node structure is where almost all of the
 *   memory allocated by this program goes.  Thus, it is imperative we
 *   keep it as small as possible so we can explore patterns as large
 *   and as deep as possible.
 *
 *   But first, how does this program even work?  Well, there are
 *   two major tricks that are used.
 *
 *   The first trick is to represent the 2D space `symbolically'
 *   (in the sense that a binary decision diagram is a symbolic
 *   representation of a boolean predicate).  This can be thought
 *   of as a sort of compression.  We break up space into a grid of
 *   squares, each containing 8x8 cells.  And we `canonicalize'
 *   each square; that is, all the squares with no cells set are
 *   represented by a single actual instance of an empty square;
 *   all squares with only the upper-left-most cell set are
 *   represented by yet another instance, and so on.  A single pointer
 *   to the single instance of each square takes less space than
 *   representing the actual cell bits themselves.
 *
 *   Where do we store these pointers?  At first, one might envision
 *   a large two-dimensional array of pointers, each one pointing
 *   to one of the square instances.  But instead, we group the
 *   squares (we'll call them 8-squares) into larger squares 16
 *   cells on a side; these are 16-squares.  Each 16-square contains
 *   four 8-squares, so each 16-square is represented by four
 *   pointers, each to an 8-square.  And we canonicalize these as
 *   well, so for a particular set of values for a 16 by 16 array
 *   of cells, we'll only have a single 16-square.
 *
 *   And so on up; we canonicalize 32-squares out of 16-squares, and
 *   on up to some limit.  Now the limit need not be very large;
 *   having just 20 levels of nodes gives us a universe that is
 *   4 * 2**20 or about 4M cells on a side.  Having 100 levels of
 *   nodes (easily within the limits of this program) gives us a
 *   universe that is 4 * 2**100 or about 5E30 cells on a side.
 *   I've run universes that expand well beyond 1E50 on a side with
 *   this program.
 *
 *   [A nice thing about this representation is that there are no
 *   coordinate values anywhere, which means that there are no
 *   limits to the coordinate values or complex multidimensional
 *   arithmetic needed.]
 *
 *   [Note that this structure so far is very similar to the octtrees
 *   used in 3D simulation and rendering programs.  It's different,
 *   however, in that we canonicalize the nodes, and also, of course,
 *   in that it is 2D rather than 3D.]
 *
 *   I mentioned there were two tricks, and that's only the first.
 *   The second trick is to cache the `results' of the LIFE calculation,
 *   but in a way that looks ahead farther in time as you go higher
 *   in the tree, much like the tree nodes themselves scan larger
 *   distances in space.  This trick is just a little subtle, but it
 *   is where the bulk of the power of the program comes from.
 *
 *   Consider once again the 8-squares.  We want to cache the result
 *   of executing LIFE on that area.  We could cache the result of
 *   looking ahead just one generation; that would yield a 6x6 square.
 *   (Note that we cannot calculate an 8-square, because we are
 *   using the single instance of the 8-square to represent all the
 *   different places that 8x8 arrangement occurs, and those different
 *   places might be surrounded by different border cells.  But we
 *   can say for sure that the central 6-square will evolve in a
 *   unique way in the next generation.)
 *
 *   We could also calculate the 4-square that is two generations
 *   hence, and the 3-square that is three generations hence, and
 *   the 2-square that is four generations hence.  We choose the
 *   4-square that is two generations hence; why will be clear in
 *   a moment.
 *
 *   Now let's consider the 16-square.  We would like to look farther
 *   ahead for this square (if we always only looked two generations
 *   ahead, our runtime would be at *least* linear in the number of
 *   generations, and we want to beat that.)  So let's look 4 generations
 *   ahead, and cache the resulting 8-square.  So we do.
 *
 *   Where do we cache the results?  Well, we cache the results in the
 *   same node structure we are using to store the pointers to the
 *   smaller squares themselves.  And since we're hashing them all
 *   together, we want a next pointer for the hash chain.  Put all of
 *   this together, and you get the following structure for the 16-squares
 *   and larger:
 */
typedef uint32_t noderef_t;
/* A noderef_t is 1 bit of GC status, 10 bits worth of offset, and 21 bits worth of allocation number. */
struct node {
   noderef_t next ;              /* hash link */
   noderef_t nw, ne, sw, se; /* constant; nw != 0 means nonleaf */
   union {
     noderef_t res ;               /* cache */
     int *resp; /* Sigh. */
   };
} ;
/*
 *   For the 8-squares, we do not have `children', we have actual data
 *   values.  We still break up the 8-square into 4-squares, but the
 *   4-squares only have 16 cells in them, so we represent them directly
 *   by an unsigned short (in this case, the direct value itself takes
 *   less memory than the pointer we might replace it with).
 *
 *   One minor trick about the following structure.  We did lie above
 *   somewhat; sometimes the struct node * points to an actual struct
 *   node, and sometimes it points to a struct leaf.  So we need a way
 *   to tell if the thing we are pointing at is a node or a leaf.  We
 *   could add another bit to the node structure, but this would grow
 *   it, and we want it to stay as small as possible.  Now, notice
 *   that, in all valid struct nodes, all four pointers (nw, ne, sw,
 *   and se) must contain a live non-zero value.  We simply ensure
 *   that the struct leaf contains a zero where the first (nw) pointer
 *   field would be in a struct node.
 *
 *   Each short represents a 4-square in normal, left-to-right then top-down
 *   order from the most significant bit.  So bit 0x8000 is the upper
 *   left (or northwest) bit, and bit 0x1000 is the upper right bit, and
 *   so on.
 */
struct leaf {
   noderef_t next ;              /* hash link */
   noderef_t isnode ;            /* must always be zero for leaves */
   unsigned short nw, ne, sw, se ;  /* constant */
   unsigned short res1, res2 ;      /* constant */
} ;
/*
 *   If it is a struct node, this returns a non-zero value, otherwise it
 *   returns a zero value.
 */
#define is_node(n) (((struct leaf *)(n))->isnode)

/*
 *   A lot of the routines from here on down traverse the universe, hanging
 *   information off the nodes.  The way they generally do so is by using
 *   (or abusing) the cache (res) field, and the least significant bit of
 *   the hash next field (as a visited bit).
 */
#define marked(n) (1 & (n)->next)
#define mark(n) ((n)->next |= 1)
#define clearmark(n) ((n)->next &= ~1)
#define clearmarkbit(p) (~1 & p)


/* The allocations come in generations -- to make it easier to look up a certain node or leaf, we put them in a union... */
union nodeleaf {
  struct leaf l;
  struct node n;
};
union nodeleaf **allocs = NULL;
int nallocs = 0;
#define ALLOCSZ 0x400
#define deref(nr)   ({noderef_t __NR = (nr); &(allocs[(__NR) >> 11][((__NR) >> 1) & 0x3FF].n);})
#define deref_l(nr) ({noderef_t __NR = (nr); &(allocs[(__NR) >> 11][((__NR) >> 1) & 0x3FF].l);})
#define ref(a, p) (((a) << 11) | ((p) << 1))

/*
 *   A key datastructure is our hash table; we use a simple bucket hash.
 */
unsigned int hashpop, hashlimit, hashprime = 100000 ;
noderef_t *hashtab ;
/*
 *   Prime hash sizes tend to work best.
 */
int nextprime(int i) {
   int j ;
   i |= 1 ;
   for (;; i+=2) {
      for (j=3; j*j<=i; j+=2)
         if (i % j == 0)
            break ;
      if (j*j > i)
         return i ;
   }
}
/*
 *   Now let's focus on some low-level details and bit manipulation.  We
 *   need a small array to hold the number of bits set in a short.  We
 *   need another small array to hold the 2-square result of a single
 *   life generation on a 4-square.  The 2-square result of that 4-square 
 *   is stored in a special way to make it easier to combine four of
 *   these into a new 4-square.  Specifically, the nw bit is stored at
 *   0x20, the ne bit at 0x10, the sw bit at 0x2, and the se bit at 0x1.
 *
 *   We initialize both of these arrays in the init() subroutine (way)
 *   below.
 *
 *   Note that all the places we represent 4-squares by short, we use
 *   unsigned shorts; this is so we can directly index into these arrays.
 */
unsigned char shortpop[65536] ;
unsigned char liferules[65536] ; /* returns four bits in the format xx..yy */
/*
 *   The cached result of an 8-square is a new 4-square representing
 *   two generations into the future.  This subroutine calculates that
 *   future, assuming liferules above is calculated.  The code that it
 *   uses is similar to code you'll see again, so we explain what's
 *   going on in some detail.
 *
 *   Each time we build a leaf node, we compute the result, because it
 *   is reasonably quick.
 *
 *   The first generation result of an 8-square is a 6-square, which
 *   we represent as nine 2-squares.  The nine 2-squares are called
 *   t00 through t22, and are arranged in a matrix:
 *
 *      t00   t01   t02
 *      t10   t11   t12
 *      t20   t21   t22
 *
 *   To compute each of these, we need to extract the relevant bits
 *   from the four 4-square values n->nw, n->ne, n->sw, and n->ne.
 *   We can use these values to directly index into the liferules
 *   array.
 *
 *   Then, given the nine values, we can compute a resulting 4-square
 *   by computing four 2-square results, and combining these into a
 *   single 4-square.
 *
 *   It's a bit intricate, but it's not really overwhelming.
 */
#define combine9(t00,t01,t02,t10,t11,t12,t20,t21,t22) \
       ((t00) << 15) | ((t01) << 13) | (((t02) << 11) & 0x1000) | \
       (((t10) << 7) & 0x880) | ((t11) << 5) | (((t12) << 3) & 0x110) | \
       (((t20) >> 1) & 0x8) | ((t21) >> 3) | ((t22) >> 5)
void leafres(struct leaf *n) {
   unsigned short
   t00 = liferules[n->nw],
   t01 = liferules[((n->nw << 2) & 0xcccc) | ((n->ne >> 2) & 0x3333)],
   t02 = liferules[n->ne],
   t10 = liferules[((n->nw << 8) & 0xff00) | ((n->sw >> 8) & 0x00ff)],
   t11 = liferules[((n->nw << 10) & 0xcc00) | ((n->ne << 6) & 0x3300) |
                   ((n->sw >> 6) & 0x00cc) | ((n->se >> 10) & 0x0033)],
   t12 = liferules[((n->ne << 8) & 0xff00) | ((n->se >> 8) & 0x00ff)],
   t20 = liferules[n->sw],
   t21 = liferules[((n->sw << 2) & 0xcccc) | ((n->se >> 2) & 0x3333)],
   t22 = liferules[n->se] ;
   n->res1 = combine9(t00,t01,t02,t10,t11,t12,t20,t21,t22) ;
   n->res2 =
   (liferules[(t00 << 10) | (t01 << 8) | (t10 << 2) | t11] << 10) |
   (liferules[(t01 << 10) | (t02 << 8) | (t11 << 2) | t12] << 8) |
   (liferules[(t10 << 10) | (t11 << 8) | (t20 << 2) | t21] << 2) |
    liferules[(t11 << 10) | (t12 << 8) | (t21 << 2) | t22] ;
}
/*
 *   Let's worry about allocation of nodes later; here's the declarations
 *   we need.  Each of these allocators actually *copy* an existing node
 *   in all of its fields.
 */
noderef_t newnode();
noderef_t newleaf();
/*
 *   We do now support garbage collection, but there are some routines we
 *   call frequently to help us.
 */
noderef_t save(noderef_t) ;
void clearstack(), pop() ;
int gsp ;
unsigned int alloced, maxmem = 2000000000 ;
#define node_hash(a,b,c,d) (((int)d)+3*(((int)c)+3*(((int)b)+3*((int)a)+3)))
#define leaf_hash(a,b,c,d) ((d)+9*((c)+9*((b)+9*(a))))
/*
 *   Resize the hash.
 */
void resize() {
   int i, nhashprime = nextprime(2 * hashprime) ;
   noderef_t p;
   noderef_t *nhashtab ;
   struct timeval t1, t2;
   if (nhashprime * sizeof(noderef_t *) > maxmem - alloced) {
      fprintf(stderr, "{no memory for new hash table}");
      hashlimit = 2000000000 ;
      return ;
   }
   /*
    *   Don't let the hash table buckets take more than 4% of the
    *   memory.  If we're starting to strain memory, let the buckets
    *   fill up a bit more.
    */
   if (nhashprime > maxmem/100) {
      fprintf(stderr, "{hash table too huge}");
      nhashprime = nextprime(maxmem/100) ;
      if (nhashprime == hashprime) {
         hashlimit = 2000000000 ;
         return ;
      }
   }
   nhashtab = (noderef_t *)calloc(nhashprime, sizeof(noderef_t)) ;
   alloced += sizeof(noderef_t) * (nhashprime - hashprime) ;
   fprintf(stderr, "{%d->%d", hashprime, nhashprime) ;
   gettimeofday(&t1, NULL);
   for (i=0; i<hashprime; i++) {
      for (p=hashtab[i]; p;) {
         struct node *ps = deref(p);
         noderef_t np = ps->next;
         unsigned int h ;
         if (is_node(ps)) {
            h = node_hash(ps->nw, ps->ne, ps->sw, ps->se) ;
         } else {
            struct leaf *l = (struct leaf *)ps ;
            h = leaf_hash(l->nw, l->ne, l->sw, l->se) ;
         }
         h %= nhashprime ;
         ps->next = nhashtab[h] ;
         nhashtab[h] = p ;
         p = np ;
      }
   }
   free(hashtab) ;
   hashtab = nhashtab ;
   hashprime = nhashprime ;
   hashlimit = hashprime / 2;
   gettimeofday(&t2, NULL);
   fprintf(stderr, "@%ldus}", (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec));
}
#ifdef TRACE
unsigned int hits ;
#endif
/*
 *   These next two routines are (nearly) our only hash table access
 *   routines; we simply look up the passed in information.  If we
 *   find it in the hash table, we return it; otherwise, we build a
 *   new node and store it in the hash table, and return that.
 */
noderef_t find_node(noderef_t nw, noderef_t ne,
                    noderef_t sw, noderef_t se) {
   struct node *p ;
   noderef_t pr;
   unsigned int h = node_hash(nw,ne,sw,se) ;
   struct node *pred = 0 ;
   h = h % hashprime ;
   for (pr=hashtab[h]; pr && (p=deref(pr)); pr = p->next) { /* make sure to compare nw *first* */
      if (nw == p->nw && ne == p->ne && sw == p->sw && se == p->se) {
#ifdef TRACE
         hits++ ;
#endif
         if (pred) { /* move this one to the front */
            pred->next = p->next ;
            p->next = hashtab[h] ;
            hashtab[h] = pr ;
         }
         return save(pr) ;
      }
      pred = p ;
   }
   pr = newnode() ;
   p = deref(pr);
   p->nw = nw ;
   p->ne = ne ;
   p->sw = sw ;
   p->se = se ;
   p->res = 0 ;
#ifdef TRACE
   if (hits > 0) {
      printf("%d hits elided\n", hits) ;
      hits = 0 ;
   }
   printf("newnode(%d %08x %08x %08x %08x) => (%08x)\n", 
                                         node_depth(p), nw, ne, sw, se, p) ;
#endif
   p->next = hashtab[h] ;
   hashtab[h] = pr ;
   hashpop++ ;
   if (hashpop > hashlimit)
      resize() ;
   return save(pr);
}
noderef_t find_leaf(unsigned short nw, unsigned short ne,
                       unsigned short sw, unsigned short se) {
   struct leaf *p ;
   noderef_t pr;
   struct leaf *pred = 0 ;
   unsigned int h = leaf_hash(nw, ne, sw, se) ;
   h = h % hashprime ;
   for (pr=hashtab[h]; pr && (p=(struct leaf *)deref(pr)); pr = p->next) {
      if (nw == p->nw && ne == p->ne && sw == p->sw && se == p->se &&
          !is_node(p)) {
#ifdef TRACE
         hits++ ;
#endif
         if (pred) {
            pred->next = p->next ;
            p->next = hashtab[h] ;
            hashtab[h] = pr ;
         }
         return save(pr) ;
      }
      pred = p ;
   }
   pr = newleaf();
   p = (struct leaf *)deref(pr);
   p->nw = nw ;
   p->ne = ne ;
   p->sw = sw ;
   p->se = se ;
   leafres(p) ;
#ifdef TRACE
   if (hits > 0) {
      printf("%d hits elided\n", hits) ;
      hits = 0 ;
   }
   printf("newleaf(%04x %04x %04x %04x) => %04x %04x (%08x)\n", nw, ne, sw, se,
                                                         p->res1, p->res2, p) ;
#endif
   p->isnode = 0 ;
   p->next = hashtab[h] ;
   hashtab[h] = pr ;
   hashpop++ ;
   if (hashpop > hashlimit)
      resize() ;
   return save(pr) ;
}
/*
 *   We've shown how to calculate the result for an 8-square.  What about
 *   the bigger squares?  Well, let's assume we have the following
 *   routines that do that work the hard way.
 */
noderef_t dorecurs_leaf() ;
noderef_t dorecurs_leaf_half() ;
noderef_t dorecurs_leaf_quarter() ;
noderef_t dorecurs() ;
noderef_t dorecurs_half() ;
int ngens = 0 ; /* must be a power of two */
int node_depth(noderef_t) ;
/*
 *   The following routine does the same, but first it checks to see if
 *   the cached result is any good.  If it is, it directly returns that.
 *   Otherwise, it figures out whether to call the leaf routine or the
 *   non-leaf routine by whether two nodes down is a leaf node or not.
 *   (We'll understand why this is a bit later.)  All the sp stuff is
 *   stack pointer and garbage collection stuff.
 */
int halvesdone ;
noderef_t getres(noderef_t nr, int depth) {
   struct node *n = deref(nr);
   if (n->res == 0) {
      int sp = gsp ;
#ifdef TRACE
      printf("getres(%08x, %d)\n", n, depth) ;
#endif
      depth-- ;
      if (ngens >= depth) {
         if (is_node(deref(n->nw))) {
            n->res = dorecurs(n->nw, n->ne, n->sw, n->se, depth) ;
         } else {
            n->res = dorecurs_leaf(n->nw, n->ne, n->sw, n->se) ;
         }
      } else {
         if (halvesdone < 1000)
            halvesdone++ ;
         if (is_node(deref(n->nw))) {
            n->res = dorecurs_half(n->nw, n->ne, n->sw, n->se, depth) ;
         } else if (ngens == 0) {
            n->res = dorecurs_leaf_quarter(n->nw, n->ne, n->sw, n->se) ;
         } else {
            n->res = dorecurs_leaf_half(n->nw, n->ne, n->sw, n->se) ;
         }
      }
      pop(sp) ;
   }
   return save(n->res);
}
/*
 *   So let's say the cached way failed.  How do we do it the slow way?
 *   Recursively, of course.  For an n-square (composed of the four
 *   n/2-squares passed in, compute the n/2-square that is n/4
 *   generations ahead.
 *
 *   This routine works exactly the same as the leafres() routine, only
 *   instead of working on an 8-square, we're working on an n-square,
 *   returning an n/2-square, and we build that n/2-square by first building
 *   9 n/4-squares, use those to calculate 4 more n/4-squares, and
 *   then put these together into a new n/2-square.  Simple, eh?
 */
noderef_t dorecurs(noderef_t nr, noderef_t ner, noderef_t tr,
                   noderef_t er, int depth) {
   int sp = gsp ;
   struct node
   *n = deref(nr),
   *ne = deref(ner),
   *t = deref(tr),
   *e = deref(er);
   noderef_t
   t00 = getres(nr, depth),
   t01 = getres(find_node(n->ne, ne->nw, n->se, ne->sw), depth),
   t02 = getres(ner, depth),
   t12 = getres(find_node(ne->sw, ne->se, e->nw, e->ne), depth),
   t11 = getres(find_node(n->se, ne->sw, t->ne, e->nw), depth),
   t10 = getres(find_node(n->sw, n->se, t->nw, t->ne), depth),
   t20 = getres(tr, depth),
   t21 = getres(find_node(t->ne, e->nw, t->se, e->sw), depth),
   t22 = getres(er, depth),
   t33 = getres(find_node(t00, t01, t10, t11), depth),
   t34 = getres(find_node(t01, t02, t11, t12), depth),
   t44 = getres(find_node(t11, t12, t21, t22), depth),
   t43 = getres(find_node(t10, t11, t20, t21), depth) ;
   nr = find_node(t33, t34, t43, t44) ;
   pop(sp) ;
   return save(nr) ;
}
/*
 *   Same as above, but we only do one step instead of 2.
 */
noderef_t dorecurs_half(noderef_t nr, noderef_t ner, noderef_t tr,
                        noderef_t er, int depth) {
   int sp = gsp ;
   struct node
   *n = deref(nr),
   *ne = deref(ner),
   *t = deref(tr),
   *e = deref(er);
   noderef_t
   t00r = getres(nr, depth),
   t01r = getres(find_node(n->ne, ne->nw, n->se, ne->sw), depth),
   t10r = getres(find_node(n->sw, n->se, t->nw, t->ne), depth),
   t11r = getres(find_node(n->se, ne->sw, t->ne, e->nw), depth),
   t02r = getres(ner, depth),
   t12r = getres(find_node(ne->sw, ne->se, e->nw, e->ne), depth),
   t20r = getres(tr, depth),
   t21r = getres(find_node(t->ne, e->nw, t->se, e->sw), depth),
   t22r = getres(er, depth) ;
   struct node
   *t00 = deref(t00r),
   *t01 = deref(t01r),
   *t10 = deref(t10r),
   *t11 = deref(t11r),
   *t02 = deref(t02r),
   *t12 = deref(t12r),
   *t20 = deref(t20r),
   *t21 = deref(t21r),
   *t22 = deref(t22r);
   if (depth > 3) {
      nr = find_node(find_node(t00->se, t01->sw, t10->ne, t11->nw),
                     find_node(t01->se, t02->sw, t11->ne, t12->nw),
                     find_node(t10->se, t11->sw, t20->ne, t21->nw),
                     find_node(t11->se, t12->sw, t21->ne, t22->nw)) ;
   } else {
      nr = find_node(find_leaf(((struct leaf *)t00)->se,
                               ((struct leaf *)t01)->sw,
                               ((struct leaf *)t10)->ne,
                               ((struct leaf *)t11)->nw),
                     find_leaf(((struct leaf *)t01)->se,
                               ((struct leaf *)t02)->sw,
                               ((struct leaf *)t11)->ne,
                               ((struct leaf *)t12)->nw),
                     find_leaf(((struct leaf *)t10)->se,
                               ((struct leaf *)t11)->sw,
                               ((struct leaf *)t20)->ne,
                               ((struct leaf *)t21)->nw),
                     find_leaf(((struct leaf *)t11)->se,
                               ((struct leaf *)t12)->sw,
                               ((struct leaf *)t21)->ne,
                               ((struct leaf *)t22)->nw)) ;
   }
   pop(sp) ;
   return save(nr) ;
}
/*
 *   If the node is a 16-node, then the constituents are leaves, so we
 *   need a very similar but still somewhat different subroutine.  Since
 *   we do not (yet) garbage collect leaves, we don't need all that
 *   save/pop mumbo-jumbo.
 */
noderef_t dorecurs_leaf(noderef_t nr, noderef_t ner, noderef_t tr, noderef_t er) {
   struct leaf
   *n = deref_l(nr),
   *ne = deref_l(ner),
   *t = deref_l(tr),
   *e = deref_l(er);
   unsigned short
   t00 = n->res2,
   t01 = deref_l(find_leaf(n->ne, ne->nw, n->se, ne->sw))->res2,
   t02 = ne->res2,
   t10 = deref_l(find_leaf(n->sw, n->se, t->nw, t->ne))->res2,
   t11 = deref_l(find_leaf(n->se, ne->sw, t->ne, e->nw))->res2,
   t12 = deref_l(find_leaf(ne->sw, ne->se, e->nw, e->ne))->res2,
   t20 = t->res2,
   t21 = deref_l(find_leaf(t->ne, e->nw, t->se, e->sw))->res2,
   t22 = e->res2 ;
   return find_leaf(deref_l(find_leaf(t00, t01, t10, t11))->res2,
                    deref_l(find_leaf(t01, t02, t11, t12))->res2,
                    deref_l(find_leaf(t10, t11, t20, t21))->res2,
                    deref_l(find_leaf(t11, t12, t21, t22))->res2) ;
}
/*
 *   Same as above but we only do two generations.
 */
#define combine4(t00,t01,t10,t11) (unsigned short)\
((((t00)<<10)&0xcc00)|(((t01)<<6)&0x3300)|(((t10)>>6)&0xcc)|(((t11)>>10)&0x33))
noderef_t dorecurs_leaf_half(noderef_t nr, noderef_t ner, noderef_t tr, noderef_t er) {
   struct leaf
   *n = deref_l(nr),
   *ne = deref_l(ner),
   *t = deref_l(tr),
   *e = deref_l(er);
   unsigned short
   t00 = n->res2,
   t01 = deref_l(find_leaf(n->ne, ne->nw, n->se, ne->sw))->res2,
   t02 = ne->res2,
   t10 = deref_l(find_leaf(n->sw, n->se, t->nw, t->ne))->res2,
   t11 = deref_l(find_leaf(n->se, ne->sw, t->ne, e->nw))->res2,
   t12 = deref_l(find_leaf(ne->sw, ne->se, e->nw, e->ne))->res2,
   t20 = t->res2,
   t21 = deref_l(find_leaf(t->ne, e->nw, t->se, e->sw))->res2,
   t22 = e->res2 ;
   return find_leaf(combine4(t00, t01, t10, t11),
                    combine4(t01, t02, t11, t12),
                    combine4(t10, t11, t20, t21),
                    combine4(t11, t12, t21, t22)) ;
}
/*
 *   Same as above but we only do one generation.
 */
noderef_t dorecurs_leaf_quarter(noderef_t nr, noderef_t ner, noderef_t tr, noderef_t er) {
   struct leaf
   *n = deref_l(nr),
   *ne = deref_l(ner),
   *t = deref_l(tr),
   *e = deref_l(er);
   unsigned short
   t00 = n->res1,
   t01 = deref_l(find_leaf(n->ne, ne->nw, n->se, ne->sw))->res1,
   t02 = ne->res1,
   t10 = deref_l(find_leaf(n->sw, n->se, t->nw, t->ne))->res1,
   t11 = deref_l(find_leaf(n->se, ne->sw, t->ne, e->nw))->res1,
   t12 = deref_l(find_leaf(ne->sw, ne->se, e->nw, e->ne))->res1,
   t20 = t->res1,
   t21 = deref_l(find_leaf(t->ne, e->nw, t->se, e->sw))->res1,
   t22 = e->res1 ;
   return find_leaf(combine4(t00, t01, t10, t11),
                    combine4(t01, t02, t11, t12),
                    combine4(t10, t11, t20, t21),
                    combine4(t11, t12, t21, t22)) ;
}
/*
 *   And that is nearly it!  Now we finish up some details.  First,
 *   allocation of nodes and leaves, but in a reasonably efficient way.
 *   For garbage collection, we need this forward declaration.
 */
int die_if_gc ;
void do_gc() ;
/*
 *   We keep free nodes in a linked list for allocation, and we allocate
 *   them 1000 at a time.
 */
int curallocblk = 0;
noderef_t freenodes = 0;
int okaytogc = 0 ;           /* only true when we're running generations. */
int totalthings = 0 ;
noderef_t newnode() {
   noderef_t r;
   if (clearmarkbit(freenodes) == 0) {
      int i ;
      
      curallocblk++;
      if (curallocblk >= nallocs) {
         if (!nallocs)
            nallocs = 1;
         nallocs *= 2;
         allocs = realloc(allocs, sizeof(union nodeleaf *) * nallocs);
         allocs[0] = NULL;
      }
      allocs[curallocblk] = calloc(ALLOCSZ, sizeof(union nodeleaf));
      freenodes = curallocblk << 11;
      alloced += ALLOCSZ * sizeof(union nodeleaf) ;
      for (i = 0; i < ALLOCSZ; i++)
         allocs[curallocblk][i].n.next = (curallocblk << 11) | ((i + 1) << 1);
      allocs[curallocblk][ALLOCSZ-1].n.next = 0;
      if (alloced > maxmem)
         fprintf(stderr, "N") ;
      totalthings += ALLOCSZ ;
   }
   if (clearmarkbit(deref(freenodes)->next) == 0 &&
       alloced + ALLOCSZ * sizeof(union nodeleaf) > maxmem &&
       okaytogc) {
      if (die_if_gc) {
         fprintf(stderr, "Dying because out of memory\n") ;
         exit(10) ;
      }
      do_gc('n') ;
   }
   r = freenodes ;
   freenodes = deref(freenodes)->next ;
   return r ;
}

noderef_t newleaf() { return newnode(); }

/*
 *   Sometimes we want the new node or leaf to be automatically cleared
 *   for us.
 */
noderef_t newclearednode() {
   noderef_t r = newnode();
   memset(deref(r), 0, sizeof(struct node));
   return r;
}
noderef_t newclearedleaf() {
   noderef_t r = newleaf();
   memset(deref(r), 0, sizeof(struct leaf));
   return r;
}
/*
 *   These are the rules we run.  Note that these rules in the birth
 *   and health arrays are *total* cells (including the middle) for
 *   birth and health; this is as opposed to the input which is only
 *   neighbor counts.
 */
int wolfram = -1 ; /* not used if -1; else 0..255 */
char *cliferules = "b3s23" ;
char birth[10] ;
char health[10] ;
/*
 *   Some globals representing our universe.  The root is the
 *   real root of the universe, and the depth is the depth of the
 *   tree where 2 means that root is a leaf, and 3 means that the
 *   children of root are leaves, and so on.  The center of the
 *   root is always coordinate position (0,0), so at startup the
 *   x and y coordinates range from -4..3; in general,
 *   -(2**depth)..(2**depth)-1.  The zeronodea is an
 *   array of canonical `empty-space' nodes at various depths.
 *   The ngens is an input parameter which is the second power of
 *   the number of generations to run.
 */
noderef_t root ;
int depth ;
noderef_t *zeronodea ;
int nzeros = 0 ;
/*
 *   Here we calculate 1-square results for 9-square inputs using
 *   our population array, and our health/birth arrays.
 */
void start_life_rules() {
   int i ;

   for (i=0; i<10; i++)
      birth[i] = health[i] = 0 ;
   for (i=0; i<65536; i++)
      liferules[i] = 0 ;
   if (wolfram >= 0) {
      for (i=0; i<0x1000; i = ((i | 0x888) + 1) & 0x1777)
         if ((i & 0x653) == 0)
            liferules[i] = 1 & (wolfram >> 
               (((i & 0x100) >> 6) | ((i & 0x20) >> 4) | ((i & 0x4) >> 2))) ;
   } else {
      char *p ;
      int isbirth = 0 ;
      for (p=cliferules; *p; p++)
         if (*p == 'b' || *p == 'B' || *p == '/')
            isbirth = 1 ;
         else if (*p == 's' || *p == 'S')
            isbirth = 0 ;
         else if ('0' <= *p && *p <= '9') {
            if (isbirth == 1) {
               birth[*p - '0'] = 1 ;
            } else {
               health[*p - '0' + 1] = 1 ;
            }
         } else {
            fprintf(stderr, "Bad life rules 2 (%c).\n", *p) ;
            exit(10) ;
         }
      for (i=0; i<0x1000; i = ((i | 0x888) + 1) & 0x1777)
         liferules[i] = (i & 0x20 ? health : birth)[shortpop[i]] ;
   }
/*
 *   We expand the 9-square to 1-square array to a 16-square to
 *   4-square array by combining results.
 */
   for (i=0; i<65536; i++)
      liferules[i] = (1 & liferules[i & 0x777]) +
                     ((1 & liferules[(i >> 1) & 0x777]) << 1) +
                     ((1 & liferules[(i >> 4) & 0x777]) << 4) +
                     ((1 & liferules[(i >> 5) & 0x777]) << 5) ;
}
/*
 *   Initialization is a little complex.
 */
void init() {
   int i ;
/*
 *   The population of one-bits in an integer is one more than the
 *   population of one-bits in the integer with one fewer bit set,
 *   and we can turn off a bit by anding an integer with the next
 *   lower integer.
 */
   for (i=1; i<65536; i++)
      shortpop[i] = shortpop[i & (i - 1)] + 1 ;
   start_life_rules() ;
   hashprime = nextprime(hashprime) ;
   hashlimit = hashprime ;
   hashtab = calloc(hashprime, sizeof(noderef_t)) ;
   alloced += hashprime * sizeof(noderef_t) ;
/*
 *   We initialize our universe to be a 16-square.
 */
   root = newclearednode() ;
   depth = 3 ;
}
/*
 *   This routine expands our universe by a factor of two, maintaining
 *   centering.  We use four new nodes, and *reuse* the root so this cannot
 *   be called after we've started hashing.
 */
void pushroot_1() {
   noderef_t tr;
   struct node *t, *r = deref(root) ;
   tr = newclearednode() ; t = deref(tr);
   t->se = r->nw ;
   r->nw = tr ;
   tr = newclearednode() ; t = deref(tr);
   t->sw = r->ne ;
   r->ne = tr ;
   tr = newclearednode() ; t = deref(tr);
   t->ne = r->sw ;
   r->sw = tr ;
   tr = newclearednode() ; t = deref(tr);
   t->nw = r->se ;
   r->se = tr ;
   depth++ ;
}
/*
 *   Return the depth of this node (2 is 8x8).
 */
int node_depth(noderef_t n) {
   int depth = 2 ;
   while (is_node(deref(n))) {
      depth++ ;
      n = deref(n)->nw ;
   }
   return depth ;
}
/*
 *   This routine returns the canonical clear space node at a particular
 *   depth.
 */
noderef_t zeronode(int depth) {
   while (depth >= nzeros) {
      int nnzeros = 2 * nzeros + 10 ;
      zeronodea = (noderef_t *)realloc(zeronodea,
                                       nnzeros * sizeof(noderef_t)) ;
      alloced += (nnzeros - nzeros) * sizeof(noderef_t) ;
      while (nzeros < nnzeros)
         zeronodea[nzeros++] = 0 ;
   }
   if (zeronodea[depth] == 0) {
      if (depth == 2) {
         zeronodea[depth] = find_leaf(0, 0, 0, 0) ;
      } else {
         noderef_t z = zeronode(depth-1) ;
         zeronodea[depth] = find_node(z, z, z, z) ;
      }
   }
   return zeronodea[depth] ;
}
/*
 *   Same, but with hashed nodes.
 */
noderef_t pushroot(noderef_t nr) {
   struct node *n = deref(nr);
   int depth = node_depth(nr) ;
   noderef_t z = zeronode(depth-1) ;
   return find_node(find_node(z, z, z, n->nw),
                    find_node(z, z, n->ne, z),
                    find_node(z, n->sw, z, z),
                    find_node(n->se, z, z, z)) ;
}
/*
 *   Here is our recursive routine to set a bit in our universe.  We
 *   pass in a depth, and walk the space.  Again, a lot of bit twiddling,
 *   but really not all that complicated.  We allocate new nodes and
 *   leaves on our way down.
 *
 *   Note that at this point our universe lives outside the hash table
 *   and has not been canonicalized, and that many of the pointers in
 *   the nodes can be null.  We'll path this up in due course.
 */
void setbit(noderef_t nr, int x, int y, int depth) {
   if (depth == 2) {
      struct leaf *l = (struct leaf *)deref(nr) ;
      if (x < 0)
         if (y < 0)
            l->sw |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
         else
            l->nw |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
      else
         if (y < 0)
            l->se |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
         else
            l->ne |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
   } else {
      int w = 1 << depth ;
      noderef_t *nptr ;
      struct node *n = deref(nr);
      depth-- ;
      if (x < 0) {
         if (y < 0)
            nptr = &(n->sw) ;
         else
            nptr = &(n->nw) ;
      } else {
         if (y < 0)
            nptr = &(n->se) ;
         else
            nptr = &(n->ne) ;
      }
      if (*nptr == 0) {
         if (depth == 2)
            *nptr = newclearedleaf() ;
         else
            *nptr = newclearednode() ;
      }
      setbit(*nptr, (x & (w - 1)) - (w >> 1), (y & (w - 1)) - (w >> 1), depth) ;
   }
}
/*
 *   Our nonrecurse top-level bit setting routine simply expands the
 *   universe as necessary to emcompass the passed-in coordinates, and
 *   then invokes the recursive setbit.
 */
void set(int x, int y) {
   while (x >= (1 << depth) || y >= (1 << depth) ||
          -x > (1 << depth) || -y > (1 << depth))
      pushroot_1() ;
   setbit(root, x, y, depth) ;
}
/*
 *   Canonicalize a universe by filling in the null pointers and then
 *   invoking find_node on each node.  Drops the original universe on
 *   the floor [big deal, it's probably small anyway].
 */
noderef_t hashpattern(noderef_t root, int depth) {
   noderef_t r;
   if (root == 0) {
      r = zeronode(depth) ;
   } else if (depth == 2) {
      struct leaf *n = (struct leaf *)deref(root) ;
      r = find_leaf(n->nw, n->ne, n->sw, n->se) ;
      n->next = freenodes ;
      freenodes = root ;
   } else {
      depth-- ;
      struct node *rootn = deref(root);
      r = find_node(hashpattern(rootn->nw, depth),
                    hashpattern(rootn->ne, depth),
                    hashpattern(rootn->sw, depth),
                    hashpattern(rootn->se, depth)) ;
      rootn->next = freenodes ;
      freenodes = root ;
   }
   return r ;
}
/*
 *   Read the RLE format.
 */
void readrle(char *line) {
   int n=0, x=0, y=0 ;
   char *p ;
   while (1) {
      if (line[0] == '#') {
         if (line[1] == 'r') {
            cliferules = strdup(line) ;
            cliferules += 2 ;
            while (*cliferules && *cliferules <= ' ')
               cliferules++ ;
            p = cliferules ;
            while (*p > ' ')
               p++ ;
            *p = 0 ;
            start_life_rules() ;
         }
      } else if (line[0] == 'x') { /* ignore x and y */
         for (p=line; *p && *p != 'r'; p++) ;
         if (strncmp(p, "rule", 4) == 0) {
            p += 4 ;
            while (*p && (*p <= ' ' || *p == '='))
               p++ ;
            cliferules = p ;
            while (*p > ' ')
               p++ ;
            *p = 0 ;
            cliferules = strdup(cliferules) ;
            start_life_rules() ;
         }
      } else {
         n = 0 ;
         for (p=line; *p; p++) {
            if ('0' <= *p && *p <= '9') {
               n = n * 10 + *p - '0' ;
            } else {
               if (n == 0)
                  n = 1 ;
               if (*p == 'b') {
                  x += n ;
               } else if (*p == 'o') {
                  while (n-- > 0)
                     set(x++, y) ;
               } else if (*p == '$') {
                  x = 0 ;
                  y -= n ;
               } else if (*p == '!') {
                  goto done ;
               } else if (*p > ' ') {
                  fprintf(stderr, "Saw illegal char %c\n", *p) ;
                  exit(10) ;
               }
               n = 0 ;
            }
         }
      }
      if (fgets(line, 10000, stdin) == 0)
         break ;
   }
done:
   root = hashpattern(root, depth) ;
}
/*
 *   This ugly bit of code will go undocumented.  It reads Alan Hensel's
 *   life format, either 1.05 or 1.06.
 */
void readpicmode(char *line) {
   int x=0, y=0 ;
   int leftx = x ;
   char *p ;
   for (;fgets(line, 10000, stdin);) {
      if (line[0] == '#') {
         if (line[1] == 'P') {
            sscanf(line + 2, " %d %d", &x, &y) ;
            leftx = x ;
            y = - y ;
         } else if (line[1] == 'N') {
            cliferules = "23/3" ;
            start_life_rules() ;
         } else if (line[1] == 'R') {
            cliferules = strdup(line) ;
            cliferules += 2 ;
            while (*cliferules && *cliferules <= ' ')
               cliferules++ ;
            p = cliferules ;
            while (*p > ' ')
               p++ ;
            *p = 0 ;
            start_life_rules() ;
         }
      } else if (line[0] == '-' || ('0' <= line[0] && line[0] <= '9')) {
         sscanf(line, "%d %d", &x, &y) ;
         set(x, y) ;
      } else if (line[0] == '.' || line[0] == '*') {
         for (p = line; *p; p++) {
            if (*p == '*')
               set(x, y) ;
            x++ ;
         }
         x = leftx ;
         y-- ;
      }
   }
   root = hashpattern(root, depth) ;
}
void readmacrocell(char *line) {
   int n=0, i=1, nw, ne, sw, se, r, d, indlen=0 ;
   noderef_t *ind = 0 ;
   while (fgets(line, 10000, stdin)) {
      if (i >= indlen) {
         int nlen = i + indlen + 10 ;
         ind = (noderef_t *)realloc(ind, sizeof(noderef_t) * nlen) ;
         while (indlen < nlen)
            ind[indlen++] = 0 ;
      }
      if (line[0] == '.' || line[0] == '*' || line[0] == '$') {
         int x=0, y=7 ;
         unsigned short lnw=0, lne=0, lsw=0, lse=0 ;
         char *p = 0 ;
         for (p=line; *p > ' '; p++) {
            switch(*p) {
case '*':      if (x > 7 || y < 0) {
                  fprintf(stderr, "Illegal coordinates (%d,%d)\n", x, y) ;
                  exit(10) ;
               }
               if (x < 4)
                  if (y < 4)
                     lsw |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
                  else
                     lnw |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
               else
                  if (y < 4)
                     lse |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
                  else
                     lne |= 1 << (3 - (x & 3) + 4 * (y & 3)) ;
               /* note: fall through here */
case '.':      x++ ;
               break ;
case '$':      x = 0 ;
               y-- ;
               break ;
default:       fprintf(stderr, "Saw illegal char %c\n", *p) ;
               exit(10) ;
            }
         }
         ind[i++] = find_leaf(lnw, lne, lsw, lse) ;
      } else {
         n = sscanf(line, "%d %d %d %d %d %d", &d, &nw, &ne, &sw, &se, &r) ;
         if (n == 0)
            continue ;
         if (n < 5) {
            fprintf(stderr, "Parse error; line is \"%s\"\n", line) ;
            exit(10) ;
         }
         if (d < 4) {
            fprintf(stderr, "Oops; depth of %d illegal here\n", d) ;
            exit(10) ;
         }
         ind[0] = zeronode(d-2) ; /* allow zeros to work right */
         if (nw < 0 || nw >= i || ind[nw] == 0 ||
             ne < 0 || ne >= i || ind[ne] == 0 ||
             sw < 0 || sw >= i || ind[sw] == 0 ||
             se < 0 || se >= i || ind[se] == 0) {
            fprintf(stderr, "Node out of range or undefined node in \"%s\"\n",
                             line) ;
            exit(10) ;
         }
         root = ind[i++] = find_node(ind[nw], ind[ne], ind[sw], ind[se]) ;
         depth = d - 1 ;
      }
   }
   if (ind)
      free(ind) ;
}
void readpattern() {
   char line[10000] ;
   if (fgets(line, 10000, stdin) == 0) {
      fprintf(stderr, "Couldn't read legal pattern file.\n") ;
      exit(10) ;
   }
   if (line[0] == '#' && line[1] == 'L') {
      readpicmode(line) ;
   } else if (line[0] == '[') {
      readmacrocell(line) ;
   } else {
      readrle(line) ;
   }
}
/*
 *   Pop off any levels we don't need.
 */
noderef_t popzeros(noderef_t nr) {
   int depth = node_depth(nr) ;
   while (depth > 3) {
      noderef_t z = zeronode(depth-2) ;
      struct node
      *n = deref(nr),
      *nw = deref(n->nw),
      *ne = deref(n->ne),
      *sw = deref(n->sw),
      *se = deref(n->se);
      if (nw->nw == z && nw->ne == z && nw->sw == z &&
          ne->nw == z && ne->ne == z && ne->se == z &&
          sw->nw == z && sw->sw == z && sw->se == z &&
          se->ne == z && se->sw == z && se->se == z) {
         depth-- ;
         nr = find_node(nw->se, ne->sw, sw->ne, se->nw) ;
      } else {
         break ;
      }
   }
   return nr;
}
/*
 *   At this point, we're pretty much done.  The remaining routines just
 *   implement various utility functions and the like, but feel free to
 *   skip down to main if you like.
 *
 *   The next few routines are a very simple multiprecision arithmetic
 *   setup using arrays of integers.
 */
int base[200] ; /* for first 100 numbers 0..99 */
/*
 *   Give me a new small multiprecision number for 0..99.
 */
int *newsmallint(int n) {
   int *r = base + 2 * n ;
   r[0] = n ;
   r[1] = -1 ;
   return r ;
}
/*
 *   Memory allocation for numbers.
 */
#define NUMMEM (16000)
struct numalloc {
   struct numalloc *next ;
   int nums[NUMMEM] ;
} *firstnummem, *currentnummem, *lastnummem ;
int numleft ;
/*
 *   Give me a pointer to a number array big enough.
 */
int *getnummem(int n) {
   if (n > NUMMEM) {
      fprintf(stderr, "Increase size of bignums") ;
      exit(10) ;
   }
   while (currentnummem == 0 || numleft < n) {
      if (currentnummem != 0 && numleft < n) {
         currentnummem = currentnummem->next ;
         numleft = NUMMEM ;
      }
      if (currentnummem == 0) {
         currentnummem = calloc(1, sizeof(struct numalloc)) ;
         alloced += sizeof(struct numalloc) ;
         if (lastnummem == 0)
            firstnummem = lastnummem = currentnummem ;
         else {
            lastnummem->next = currentnummem ;
            lastnummem = currentnummem ;
         }
         numleft = NUMMEM ;
      }
   }
   numleft -= n ;
   return currentnummem->nums + numleft ;
}
/*
 *   Reset num memories.
 */
void resetnummem() {
   currentnummem = firstnummem ;
}
/*
 *   Given an accumulator which is big enough, add s.
 */
void checkgensize(int) ;
int gensize ;
int addto(int *acc, int *s) {
   int r = 0 ;
   int m = 1 ;
   while (*acc != -1 || *s != -1 || r > 0) {
      if (*acc == -1) {
         *acc = 0 ;
         acc[1] = -1 ;
      }
      if (*s == -1)
         *acc += r ;
      else
         *acc += *s++ + r ;
      if (*acc >= 1000000000) {
         r = 1 ;
         *acc -= 1000000000 ;
      } else
         r = 0 ;
      acc++ ;
      m++ ;
   }
   if (m >= gensize)
      checkgensize(m) ;
   return m ;
}
/*
 *   Return the sum of four numbers.
 */
int *tsum4 ;
int *sum4(int *a1, int *a2, int *a3, int *a4) {
   int m = 0 ;
   int i ;
   int *r ;
   for (i=0; a1[i] != -1; i++)
      tsum4[i] = a1[i] ;
   tsum4[i] = -1 ;
   addto(tsum4, a2) ;
   addto(tsum4, a3) ;
   m = addto(tsum4, a4) ;
   if (m == 2 && tsum4[0] < 100)
      return newsmallint(tsum4[0]) ;
   else {
      r = getnummem(m) ;
      memcpy(r, tsum4, m * sizeof(int)) ;
      return r ;
   }
}
/*
 *   Stringify one of these multiprecision numbers.
 */
char *tstringify ;
char *stringify(int *n) {
   int i ;
   char *r = tstringify ;
   for (i=0; n[i+1] != -1; i++) ;
   if (n[i] > 999999) {
      sprintf(r, "%d,%03d,%03d", n[i]/1000000, (n[i]/1000)%1000, n[i]%1000) ;
   } else if (n[i] > 999) {
      sprintf(r, "%d,%03d", n[i]/1000, n[i]%1000) ;
   } else {
      sprintf(r, "%d", n[i]) ;
   }
   while (--i >= 0) {
      r += strlen(r) ;
      sprintf(r, ",%03d,%03d,%03d", n[i]/1000000, (n[i]/1000)%1000, n[i]%1000) ;
   }
   return tstringify ;
}
/*
 *   This recursive routine calculates the population by hanging the
 *   population on marked nodes.
 */
int *calcpop(noderef_t root, int depth) {
   int *r ;
   if (root == zeronode(depth)) {
      return newsmallint(0) ;
   } else if (depth == 2) {
      struct leaf *n = (struct leaf *)deref(root) ;
      return newsmallint(shortpop[n->nw] + shortpop[n->ne] +
                         shortpop[n->sw] + shortpop[n->se]) ;
   } else if (marked(deref(root))) {
      return deref(root)->resp ;
   } else {
      struct node *rr = deref(root);
      depth-- ;
      r = sum4(calcpop(rr->nw, depth), calcpop(rr->ne, depth),
               calcpop(rr->sw, depth), calcpop(rr->se, depth)) ;
      mark(deref(root)) ;
      rr->resp = r ;
      return r ;
   }
}
/*
 *   After running one of the marking routines, we can clean up by
 *   calling this routine.
 */
void aftercalcpop(noderef_t rootr, int depth) {
   struct node *root = deref(rootr);
   if (marked(root)) {
      if (depth == 2) {
         root->nw = 0 ;
         clearmark(root) ;
      } else {
         root->resp = 0 ;  /* clear cache field! */
         clearmark(root) ;
         depth-- ;
         aftercalcpop(root->nw, depth) ;
         aftercalcpop(root->ne, depth) ;
         aftercalcpop(root->sw, depth) ;
         aftercalcpop(root->se, depth) ;
      }
   }
}
/*
 *   This top level routine calculates the population of a universe.
 */
char *population(noderef_t root) {
   int depth, *r ;
   resetnummem() ;
   depth = node_depth(root) ;
   r = calcpop(root, depth) ;
   aftercalcpop(root, depth) ;
   return stringify(r) ;
}
/*
 *   Unpack into two 8x8s into a more conventional arrangement.
 */
void unpack8x8(unsigned short nw, unsigned short ne, unsigned short sw,
               unsigned short se, unsigned int *top, unsigned int *bot) {
   *top = ((nw & 0xf000) << 16) | (((ne & 0xf000) | (nw & 0xf00)) << 12) |
          (((ne & 0xf00) | (nw & 0xf0)) << 8) |
          (((ne & 0xf0) | (nw & 0xf)) << 4) | (ne & 0xf) ;
   *bot = ((sw & 0xf000) << 16) | (((se & 0xf000) | (sw & 0xf00)) << 12) |
          (((se & 0xf00) | (sw & 0xf0)) << 8) |
          (((se & 0xf0) | (sw & 0xf)) << 4) | (se & 0xf) ;
}
/*
 *   This routine recursively writes out the cell structure of the universe.
 */
int cellcounter = 0 ;
FILE *macro ;
void writecell(noderef_t rootr, int depth) {
   int thiscell = 0 ;
   struct node *root = deref(rootr);
   if (marked(root))
      return ;
   mark(root) ;
   if (rootr == zeronode(depth)) {
      if (depth == 2)
         root->nw = 0 ;
      else
         root->res = 0 ;
   } else if (depth == 2) {
      int i, j ;
      unsigned int top, bot ;
      struct leaf *n = (struct leaf *)root ;
      thiscell = ++cellcounter ;
      root->nw = (noderef_t) thiscell ;
      unpack8x8(n->nw, n->ne, n->sw, n->se, &top, &bot) ;
      for (j=7; (top | bot) && j>=0; j--) {
         int bits = (top >> 24) ;
         top = (top << 8) | (bot >> 24) ;
         bot = (bot << 8) ;
         for (i=0; bits && i<8; i++, bits = (bits << 1) & 255)
            if (bits & 128)
               fprintf(macro, "*") ;
            else
               fprintf(macro, ".") ;
         fprintf(macro, "$") ;
      }
      fprintf(macro, "\n") ;
   } else {
      writecell(root->nw, depth-1) ;
      writecell(root->ne, depth-1) ;
      writecell(root->sw, depth-1) ;
      writecell(root->se, depth-1) ;
      thiscell = ++cellcounter ;
      root->res = (noderef_t) thiscell ;
      if (depth > 3) {
         fprintf(macro, "%d %d %d %d %d\n", depth+1,
                         (int)deref(root->nw)->res, (int)deref(root->ne)->res,
                         (int)deref(root->sw)->res, (int)deref(root->se)->res) ;
      } else
         fprintf(macro, "%d %d %d %d %d\n", depth+1,
                         (int)deref(root->nw)->nw, (int)deref(root->ne)->nw,
                         (int)deref(root->sw)->nw, (int)deref(root->se)->nw) ;
   }
}
int *maxgens, *incgens, *gen, *tmain ;
/*
 *   And this is the top-level writer of the cell structure.
 */
void writeout(char *filename, int argc, char *argv[]) {
   int depth = node_depth(root) ;
   int i ;
   macro = fopen(filename, "w") ;
   if (!macro) {
      fprintf(stderr, "Couldn't open output file %s\n", filename) ;
      exit(10) ;
   }
   fprintf(stderr, "[->%s", filename) ;
   fprintf(macro, "[M1] (hlife 0.96)") ;
   for (i=1; i<argc; i++)
      fprintf(macro, " %s", argv[i]) ;
   fprintf(macro, "\n") ;
   cellcounter = 0 ;
   writecell(root, depth) ;
   fclose(macro) ;
   aftercalcpop(root, depth) ;
   fprintf(stderr, "]") ;
}
/*
 *   Finally, our gc routine.  We keep a `stack' of all the `roots'
 *   we want to preserve.  Nodes not reachable from here, we allow to
 *   be freed.  Same with leaves.
 */
noderef_t *stack ;
int stacksize ;
/*
 *   This routine marks a node as needed to be saved.
 */
noderef_t save(noderef_t n) {
   if (gsp >= stacksize) {
      int nstacksize = stacksize * 2 + 100 ;
      alloced += sizeof(noderef_t)*(nstacksize-stacksize) ;
      stack = realloc(stack, nstacksize * sizeof(noderef_t)) ;
      stacksize = nstacksize ;
   }
   stack[gsp++] = n ;
   return n ;
}
/*
 *   This routine pops the stack back to a previous depth.
 */
void pop(int n) {
   gsp = n ;
}
/*
 *   This routine clears the stack altogether.
 */
void clearstack() {
   gsp = 0 ;
}
/*
 *   Do a gc.  Walk down from all nodes reachable on the stack, saveing
 *   them by setting the odd bit on the next link.  Then, walk the hash,
 *   eliminating the res from everything that's not saveed, and moving
 *   the nodes from the hash to the freelist as appropriate.  Finally,
 *   walk the hash again, clearing the low order bits in the next pointers.
 */
void gc_mark(noderef_t rootr) {
   struct node *root = deref(rootr);
   if (!marked(root)) {
      mark(root) ;
      if (is_node(root)) {
         gc_mark(root->nw) ;
         gc_mark(root->ne) ;
         gc_mark(root->sw) ;
         gc_mark(root->se) ;
         if (root->res)
            gc_mark(root->res) ;
      }
   }
}
void do_gc(int why) {
   int i, j;
   unsigned int freed_nodes=0;
   union nodeleaf *n;
   struct timeval t1, t2;
   gettimeofday(&t1, NULL);
   fprintf(stderr, "[%c%d", why, hashpop/(totalthings/100)) ;
   for (i = 1; i <= curallocblk; i++)
      for (j = 0; j < ALLOCSZ; j++)
         clearmark(&(allocs[i][j].l));
   for (i=0; i<gsp; i++)
      gc_mark(stack[i]) ;
   hashpop = 0 ;
   memset(hashtab, 0, sizeof(noderef_t) * hashprime) ;
   fprintf(stderr, ":") ;
   freenodes = 0 ;
   for (i = 1; i <= curallocblk; i++)
      for (j = 0; j < ALLOCSZ; j++) {
         if (marked(&(allocs[i][j].l))) {
            int h;
            if (allocs[i][j].l.isnode)
              h = leaf_hash(allocs[i][j].l.nw, allocs[i][j].l.ne, allocs[i][j].l.sw, allocs[i][j].l.se) % hashprime ;
            else
              h = node_hash(allocs[i][j].n.nw, allocs[i][j].n.ne, allocs[i][j].n.sw, allocs[i][j].n.se) % hashprime ;
            allocs[i][j].n.next = hashtab[h] ;
            hashtab[h] = ref(i, j) ;
            hashpop++ ;
         } else {
            allocs[i][j].n.next = freenodes ;
            freenodes = ref(i, j) ;
            freed_nodes++ ;
         }
         clearmark(&(allocs[i][j].l));
      }
   gettimeofday(&t2, NULL);
   fprintf(stderr, "%u@%ldus]", hashpop/(totalthings/100), (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec));
}
/*
 *   Clear the cache bits down to the appropriate level, marking the
 *   nodes we've handled.
 */
void clearcache(noderef_t nr, int depth, int clearto) {
   struct node *n = deref(nr);
   if (!marked(n)) {
      mark(n) ;
      if (depth > 3) {
         depth-- ;
         clearcache(n->nw, depth, clearto) ;
         clearcache(n->ne, depth, clearto) ;
         clearcache(n->sw, depth, clearto) ;
         clearcache(n->se, depth, clearto) ;
         if (n->res)
            clearcache(n->res, depth, clearto) ;
      }
      if (depth >= clearto)
         n->res = 0 ;
   }
}
/*
 *   Change the ngens value.  Requires us to walk the hash, clearing
 *   the cache fields of any nodes that do not have the appropriate
 *   values.
 */
void new_ngens(int newval) {
   int i, n, j ;
   noderef_t p, pp;
   struct node *sp;
   int clearto = ngens ;
   if (newval > ngens && halvesdone == 0) {
      ngens = newval ;
      return ;
   }
   fprintf(stderr, "<") ;
   if (newval < clearto)
      clearto = newval ;
   clearto++ ; /* clear this depth and above */
   if (clearto < 3)
      clearto = 3 ;
   ngens = newval ;
   for (i=0; i<hashprime; i++)
      for (p=hashtab[i]; p && (sp = deref(p)); p=clearmarkbit(sp->next))
         if (is_node(sp) && !marked(sp))
            clearcache(p, node_depth(p), clearto) ;
   fprintf(stderr, "%d", halvesdone) ;
   for (i = 1; i <= curallocblk; i++)
      for (j = 0; j < ALLOCSZ; j++)
         clearmark(&(allocs[i][j].n));
   fprintf(stderr, ">") ;
   halvesdone = 0 ;
}
/*
 *   Return log2.
 */
int log2(unsigned int n) {
   int r = 0 ;
   while ((n & 1) == 0) {
      n >>= 1 ;
      r++ ;
   }
   if (n != 1) {
      fprintf(stderr, "Expected power of two!") ;
      exit(10) ;
   }
   return r ;
}
/*
 *   Check to see if we need to increase the size of our multiprecision
 *   constants.  (These are not allocated out of nummem because we
 *   reclaim all of nummem after every big calculation.)
 */
void biggergennums() {
   int ngensize = 2 * gensize + 5 ;
   maxgens = (int *)realloc(maxgens, ngensize * sizeof(int)) ;
   incgens = (int *)realloc(incgens, ngensize * sizeof(int)) ;
   gen = (int *)realloc(gen, ngensize * sizeof(int)) ;
   tsum4 = (int *)realloc(tsum4, ngensize * sizeof(int)) ;
   tmain = (int *)realloc(tmain, ngensize * sizeof(int)) ;
   tstringify = (char *)realloc(tstringify, ngensize * 12) ;
   gensize = ngensize - 2 ;
}
/*
 *   Increase gen size as needed.
 */
void checkgensize(int n) {
   while (n >= gensize)
      biggergennums() ;
}
/*
 *   Compare two big numbers.
 */
int bignumcmp(int *a, int *b) {
   int i ;
   for (i=0; a[i] != -1 && b[i] != -1; i++) ;
   if (a[i] != -1)
      return 1 ;
   if (b[i] != -1)
      return -1 ;
   for (i--; a[i] == b[i] && i > 0; i--) ;
   if (a[i] > b[i])
      return 1 ;
   if (a[i] < b[i])
      return -1 ;
   return 0 ;
}
/*
 *   Find the log base 2 of a bignum.
 */
int bignumlog2(int *p) {
   int r = 0 ;
   tmain[0] = 1 ;
   tmain[1] = -1 ;
   while (1) {
      int c = bignumcmp(p, tmain) ;
      if (c == 0)
         return r ;
      if (c < 0) {
         fprintf(stderr, "Must give power of two; saw %s\n", stringify(p)) ;
         exit(10) ;
      }
      addto(tmain, tmain) ;
      r++ ;
   }
}
/*
 *   Divide a bignum by 2.
 */
void bignumdiv2(int *p) {
   int r=0, i ;
   for (i=0; p[i] != -1; i++) ;
   i-- ;
   while (i >= 0) {
      int n = r * 500000000 + p[i] / 2 ;
      r = p[i] & 1 ;
      p[i--] = n ;
   }
   for (i=0; p[i] != -1; i++) ;
   i-- ;
   if (p[i] == 0 && i > 0)
      p[i] = -1 ;
}
/*
 *   Parse an integer.  Either decimal notation or 2^ notation.
 */
int *parseint(int **what, char *s) {
   int i, *p, v ;
   if (*s == '2' && s[1] == '^') {
      if (sscanf(s+2, "%d", &i) != 1)
         return 0 ;
      checkgensize(i/30 + 1) ;
      p = *what ;
      p[0] = 1 ;
      p[1] = -1 ;
      while (i > 0) {
         addto(*what, *what) ;
         i-- ;
      }
      return *what ;
   }
   i = strlen(s) ;
   while (gensize * 9 <= i + 9)
      biggergennums() ;
   p = *what ;
   v = 0 ;
   if (i > 1 && *s == '0')
      return 0 ;
   p[(i + 8) / 9] = -1 ;
   while (i > 0) {
      if ('0' <= *s && *s <= '9')
         v = 10 * v + *s++ - '0' ;
      else
         return 0 ;
      i-- ;
      if (i % 9 == 0) {
         p[i / 9] = v ;
         v = 0 ;
      }
   }
   return p ;
}
/*
 *   Finally, we get to run the pattern.  We first ensure that all
 *   clearspace nodes and the input pattern is never garbage
 *   collected; we turn on garbage collection, and then we invoke our
 *   magic top-level routine passing in clearspace borders that are
 *   guaranteed large enough.
 */
noderef_t runpattern(noderef_t n) {
   int depth = node_depth(n) ;
   noderef_t n2 ;
   n = pushroot(n) ;
   depth++ ;
   n = pushroot(n) ;
   depth++ ;
   while (ngens + 2 > depth) {
      n = pushroot(n) ;
      depth++ ;
   }
   okaytogc = 1 ;
   save(zeronode(nzeros-1)) ;
   save(n) ;
   n2 = getres(n, depth) ;
   okaytogc = 0 ;
   clearstack() ;
   if (halvesdone == 1) {
      deref(n)->resp = 0 ;
      halvesdone = 0 ;
   }
   n = popzeros(n2) ;
   addto(gen, incgens) ;
   return n ;
}
/*
 *   Finally, our main routine!
 */
int main(int argc, char *argv[]) {
   char *outfile = 0 ;
   int nopop = 0 ;
   int opts = 0 ; /* what options have been specified? */
   int oargc = argc ;
   char **oargv = argv ;
   printf("This is hlife 0.96 Copyright 2001-2005 Radical Eye Software\n") ;
   biggergennums() ;
   maxgens[0] = 0 ;
   maxgens[1] = -1 ;
   incgens[0] = 1 ;
   incgens[1] = -1 ;
   gen[0] = 0 ;
   gen[1] = -1 ;
   init() ;
   while (argc > 1 && argv[1][0] == '-') {
      argc-- ;
      argv++ ;
      switch (argv[0][1]) {
case '2':
         opts |= 1 << ('t' - 'a') ;
         break ;
/*
 *   Currently -M only indicates the max mem you *want* it to use, and
 *   we calculate the size of our hash table based on that.
 */
case 'M':
         if (sscanf(argv[1], "%u", &maxmem) != 1) {
            fprintf(stderr, "Need valid integer for -M (maxmem); saw %s\n",
                            argv[1]) ;
            exit(10) ;
         }
         maxmem *= 1000000 ;
         argv++ ;
         argc-- ;
         break ;
/*
 *   This argument is the interval to use.
 *   Currently must be a power of two.  Can be expressed as 2^n.
 */
case 'i':
         if (parseint(&incgens, argv[1]) == 0) {
            fprintf(stderr, "Need valid integer for -i (interval); saw %s\n",
                            argv[1]) ;
            exit(10) ;
         }
         argv++ ;
         argc-- ;
         opts |= 1 << ('i' - 'a') ;
         break ;
/*
 *   This argument is the maximum generation to compute.  Can be
 *   expressed as 2^n.
 */
case 'm':
         if (parseint(&maxgens, argv[1]) == 0) {
            fprintf(stderr, "Need valid integer for -m (maxgens); saw %s\n",
                            argv[1]) ;
            exit(10) ;
         }
         argv++ ;
         argc-- ;
         opts |= 1 << ('m' - 'a') ;
         break ;
/*
 *   This is the name of the output file to write a macrocell to.
 */
case 'o':
         outfile = argv[1] ;
         argv++ ;
         argc-- ;
         break ;
/*
 *   No population counts.
 */
case 'q':
         nopop++ ;
         break ;
/*
 *   Wolfram values.
 */
case 'w':
         if (sscanf(argv[1], "%d", &wolfram) != 1) {
            fprintf(stderr, "Need valid integer for -w (wolfram); saw %s\n",
                            argv[1]) ;
            exit(10) ;
         }
         argv++ ;
         argc-- ;
         break ;
/*
 *   Life/death rules.
 */
case 'r':
         cliferules = argv[1] ;
         start_life_rules() ;
         argv++ ;
         argc-- ;
         break ;
/*
 *   Disable gc; die if we run out of RAM.
 */
case 'f':
         die_if_gc = 1 ;
         break ;
default:
         fprintf(stderr, "Didn't understand arg %s\n", argv[0]) ;
         exit(10) ;
         break ;
      }
   }
   if (argc > 1 && freopen(argv[1], "r", stdin) == NULL) {
      fprintf(stderr, "Couldn't open %s for input\n", argv[1]) ;
      exit(10) ;
   }
/*
 *   Now we initialize things, read the pattern, clear out some space at
 *   the bottom and left (for the universe to expand; run pattern takes
 *   care of the top and right).  We then make the universe big enough for
 *   the number of generations we run (this also implicitly sets the
 *   number of generations to run because currently we basically just
 *   advance the top node), and then we run the pattern.  Finally, we
 *   calculate and print the population, and fulfill any other requests
 *   that might have been made of us.
 */
   readpattern() ;
   if (!nopop)
      printf("Start population is %s\n", population(root)) ;
   if (opts == 0 || opts & (1 << ('t' - 'a'))) { /* do 1, 2, 4, 8, 16, 32, 64, 128, ... */
#ifdef UI
      UI(root) ;
      return 0 ;
#else
      int fc = -1 ;
      char outname[100] ;
      if (!(opts & (1 << ('i' - 'a')))) {
         incgens[0] = 1 ;
         incgens[1] = -1 ;
      }
      ngens = bignumlog2(incgens) ;
      while (1) {
         root = runpattern(root) ;
         printf(">> %s", stringify(gen)) ;
         if (nopop)
            printf("\n") ;
         else
            printf(" %s\n", population(root)) ;
         fflush(stdout) ;
         fc++ ;
         if (outfile) {
            sprintf(outname, "%s-%d", outfile, fc) ;
            writeout(outname, oargc, oargv) ;
         }
         if (fc > 0) {
            new_ngens(ngens + 1) ;
            addto(incgens, incgens) ;
         }
         if ((opts & (1 << ('m' - 'a'))) &&
             bignumcmp(gen, maxgens) >= 0)
            break ;
      }
#endif
   } else if (opts == (1 << ('m' - 'a'))) { /* go to max. */
      ngens = 0 ;
      memcpy(tmain, maxgens, sizeof(int)*(gensize+2)) ;
      incgens[0] = 1 ;
      incgens[1] = -1 ;
      while (bignumcmp(tmain, newsmallint(0)) > 0) {
         int nngens = ngens ;
         while ((tmain[0] & 1) == 0) {
            bignumdiv2(tmain) ;
            addto(incgens, incgens) ;
            nngens++ ;
         }
         new_ngens(nngens) ;
         root = runpattern(root) ;
         tmain[0]-- ;
         printf(">> %s", stringify(gen)) ;
         if (nopop)
            printf("\n") ;
         else
            printf(" %s\n", population(root)) ;
         fflush(stdout) ;
      }
   } else if (opts & (1 << ('i' - 'a'))) { /* do by inc; must be power of 2 */
      ngens = bignumlog2(incgens) ;
      while (1) {
         root = runpattern(root) ;
         printf(">> %s", stringify(gen)) ;
         if (nopop)
            printf("\n") ;
         else
            printf(" %s\n", population(root)) ;
         fflush(stdout) ;
         if ((opts & (1 << ('m' - 'a'))) &&
             bignumcmp(gen, maxgens) >= 0)
            break ;
      }
   }
   if (outfile)
      writeout(outfile, oargc, oargv) ;
#ifdef UI
   UI(root) ;
#endif
/* printf("Hash pop %d\n", hashpop) ; */
   return 0 ;
}
