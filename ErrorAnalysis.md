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