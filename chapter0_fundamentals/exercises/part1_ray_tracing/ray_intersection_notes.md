# Ray–Segment Intersection: The Formula, No-Solution Cases, and the Determinant

## The formula

We want a point that is **on the ray** and **on the segment** at the same time:

$$
O + uD = L_1 + v(L_2 - L_1)
$$

Move the unknowns ($u$, $v$) to the left and the known points to the right:

$$
uD + v(L_1 - L_2) = L_1 - O
$$

Split into the $x$ and $y$ components:

$$
\begin{aligned}
u\,D_x + v\,(L_1 - L_2)_x &= (L_1 - O)_x \\
u\,D_y + v\,(L_1 - L_2)_y &= (L_1 - O)_y
\end{aligned}
$$

In matrix form, $A\mathbf{x} = \mathbf{b}$:

$$
A = \begin{pmatrix} D_x & (L_1 - L_2)_x \\ D_y & (L_1 - L_2)_y \end{pmatrix}, \qquad
\mathbf{x} = \begin{pmatrix} u \\ v \end{pmatrix}, \qquad
\mathbf{b} = \begin{pmatrix} (L_1 - O)_x \\ (L_1 - O)_y \end{pmatrix}
$$

The **columns of $A$** are the two direction vectors:

- column 1 is $D$, the direction the ray travels
- column 2 is $L_1 - L_2$, the direction the segment runs

---

## When is there no solution?

There's no unique solution when $D$ and $L_1 - L_2$ are **parallel**, meaning the ray and the segment point the same way.

Two lines in a plane can relate in three ways:

```
1. They cross          2. Parallel, apart       3. On top of each other
                                                   (collinear)
     \   /                ───────                  ═══════════
      \ /                                          (same line)
       X                  ───────
      / \
   one crossing        zero crossings           infinitely many
   → one (u, v)        → no (u, v)              → infinite (u, v)
```

Cases 2 and 3 both happen when $D \parallel (L_1 - L_2)$. In neither case is there one specific $(u, v)$ for `solve` to return.

### Example: parallel lines that never meet

$$
O = (0, 0), \quad D = (0, 1), \quad L_1 = (1, -1), \quad L_2 = (1, 1)
$$

Then $L_1 - L_2 = (0, -2)$ and $L_1 - O = (1, -1)$:

$$
\begin{pmatrix} 0 & 0 \\ 1 & -2 \end{pmatrix}\begin{pmatrix} u \\ v \end{pmatrix} = \begin{pmatrix} 1 \\ -1 \end{pmatrix}
$$

Row 1 says:

$$
0 \cdot u + 0 \cdot v = 1 \quad\Longrightarrow\quad 0 = 1
$$

That's impossible, so no $u, v$ exist. The ray stays at $x = 0$ forever and the segment sits at $x = 1$. Neither direction vector has an $x$ part, so they never meet.

### Example: on top of each other

$$
O = (0, 0), \quad D = (0, 1), \quad L_1 = (0, 1), \quad L_2 = (0, 2)
$$

Then $L_1 - L_2 = (0, -1)$ and $L_1 - O = (0, 1)$:

$$
\begin{pmatrix} 0 & 0 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} u \\ v \end{pmatrix} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}
$$

Row 1 says $0 = 0$, which is always true and gives no information. Row 2 says $u - v = 1$. That's one equation with two unknowns, so there are infinitely many solutions:

$$
(u, v) = (1, 0), \quad (1.5, 0.5), \quad (2, 1), \quad \dots
$$

ARENA says to treat **both** cases as no intersection (`False`).

---

## Where the determinant comes in

For a $2 \times 2$ matrix:

$$
\det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad - bc
$$

**Key fact: $\det(A) = 0$ exactly when the two columns are parallel.**

To see why, suppose column 2 is a scaled copy of column 1, so $(b, d) = k\,(a, c)$. Then:

$$
ad - bc = a(kc) - (ka)c = 0
$$

Check it on the examples:

$$
\det\begin{pmatrix} 0 & 0 \\ 1 & -2 \end{pmatrix} = 0 \cdot (-2) - 0 \cdot 1 = 0 \qquad \text{(parallel, apart)}
$$

$$
\det\begin{pmatrix} 0 & 0 \\ 1 & -1 \end{pmatrix} = 0 \cdot (-1) - 0 \cdot 1 = 0 \qquad \text{(collinear)}
$$

And on the earlier hit example, $D = (1, 0.5)$ with $L_1 = (1,-1)$, $L_2 = (1,1)$:

$$
\det\begin{pmatrix} 1 & 0 \\ 0.5 & -2 \end{pmatrix} = 1 \cdot (-2) - 0 \cdot 0.5 = -2 \neq 0 \qquad \text{(unique solution)}
$$

### The geometric view

$\det(A)$ is the signed **area of the parallelogram** formed by the two columns $D$ and $L_1 - L_2$:

```
  det ≠ 0                         det = 0

     L1−L2                        ────D────>
       ↗ ─────────↗               ──L1−L2──>
      /          /                (flat: no area)
     / area>0  /
    ↗─────────↗
        D
```

If the two directions differ, they span a real parallelogram with nonzero area. If they're parallel, it's flattened to a line and the area is $0$.

### Why this helps

$$
\det(A) \neq 0 \;\Longrightarrow\; \text{exactly one } (u, v) \;\Longrightarrow\; \text{solve, then check } u \ge 0 \text{ and } 0 \le v \le 1
$$

$$
\det(A) = 0 \;\Longrightarrow\; \text{parallel lines (zero or infinitely many solutions)} \;\Longrightarrow\; \text{return False}
$$

Two ways to handle the parallel case in code:

1. **`try/except`:** call `t.linalg.solve`. If $\det(A) = 0$ it raises an error, and you catch it and return `False`.
2. **Check the determinant first:** compute `t.linalg.det(A)` and return `False` if it's close to $0$, *before* calling `solve`.

For `intersect_ray_1d` either works. For the **batched** version you need option 2, because one singular matrix makes `solve` raise an error for the whole batch.

Because of floating point, a "parallel" determinant may come out as something like $10^{-9}$ instead of exactly $0$, so in practice you check

$$
\lvert \det(A) \rvert < \varepsilon
$$

for some small $\varepsilon$, rather than testing $\det(A) = 0$ exactly.
