# Error Analysis

What ends up helping the error go down? considering the strong gaussian in the center (reducing any chance of triangle interaction)
- Increasing quad refinement levels
- Increasing quadrature order

However, even with high quadrature order, the error is still limited by the center. Consult the following image:

![Order 17 domain quadrature error graph](Screenshot_20260727_234720.png)

No matter what changes are made though, the center spike doesn't really go down in error. I'm inclined to think that this is floating point error.

Additional note: the circles on the graph above are somewhat related to the near-singular potential evaluation; as I reduce the radius of target points included in the near-singular potential evaluation, these circles shrink.

# New Error Notes

Right now, I have a way to get the error down near 2e-12

8425 - 26.51 seconds - 317.8
9288 - 28.59 - 324.9
25982 - 78.86 - 329.5

Effectively, the list of knobs that I can turn right now:
- Quad Qorder
- Triangle Qorder
- Quad legendre matrix degree
- Boundary discretization points
- Boundary curvature
- Quad refinement tolerances
- Triangle refinement tolerances

to discuss:
- FMM
- better accuracy yay!
- some results
- boundary limitations + current progress
- next steps (further error reduction, tolerance function, legendre tail coeff refinement)

Reformat to use sparse correction matrix


0.535647 seconds total
0.393216 seconds FMM
0.029228 seconds N
0.039794 seconds q_map
0.128781 seconds t_map

more accurate look, through profview:
of 6614 ticks in the map application:
4508 are in the main FMM application (quads -> quads/triangles)
1721 are in the other FMM application (triangles -> quads/triangles)
311 are in the quad -> quad correction
Cannot even find the quad -> triangle application (negligible)

allocations:
fmm map: 76 allocations
N map: 3 allocations
q_map: 1.84k allocations
t_map: 140 allocations

this is pretty great!

todo still:
- test with larger num of target pts to see how it scales
- optimize a lil further, improve code qual
- check over operator creation once more lol
- test with lower error again



Accuracy vs. order

| Total cache hits | F cache hit number | F cache hit % | L cache hit number | L cache hit % |
| ---------------- | ------------------ | ------------- | ------------------ | ------------- |
| $185680$         | $183572$           | $98.8$%       | $130906$           | $70.5$%       |

note - first application usually takes longer, has more allocations. Taking other allocations here

| qorder | np  | max error | application time | application allocations |
| ------ | --- | --------- | ---------------- | ----------------------- |
| $17$   | $9$ | $4.3e-11$ | $0.52$s          | $221$                   |
| $20$   | $11$| $4.3e-11$ | $0.68$s          | $221$                   |
| $30$   | $16$| $4.3e-11$ | $1.56$s          | $222$

Increasing qorder further doesn't reduce max error but reduces L2 error, not by much tho

error went up when triangle qorder went up???
