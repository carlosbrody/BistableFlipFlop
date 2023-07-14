module BistableDynamics

export decisionMake

using PyPlot, Statistics, LinearAlgebra, Printf


###########################################
#
#   UTILITIES (AXISMOVE, ETC)
#
###########################################


"""
    ax = axisWidthChange(factor; lock="c", ax=nothing)

    Changes the width of the current axis by multiplicative factor.
    lock can be one of "c" or "m" (meaning the center of the axis
    stays fixed) or "l" (the left edge stays fixed) or "r" (the right
    edge stays fixed).

    If ax is passed then the indicated axis is modified, not necessarily
    the current one.
"""
function axisWidthChange(factor; lock="c", ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax.get_position().bounds

    if lock=="l";
    elseif lock=="c" || lock=="m"; x = x + w*(1-factor)/2;
    elseif lock=="r"; x = x + w*(1-factor);
    else error("I don't know lock type ", lock)
    end

    w = w*factor;
    ax.set_position([x, y, w, h])

    return ax
end


"""
ax = axisHeightChange(factor; lock="c", ax=nothing)

    Changes the height of the current axis by multiplicative factor.
    lock can be one of "c" or "m" (meaning the center of the axis
    stays fixed) or "b" (the bottom edge stays fixed) or "t" (the top
    edge stays fixed).

    If ax is passed then the indicated axis is modified, not necessarily
    the current one.
"""
function axisHeightChange(factor; lock="c", ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax.get_position().bounds

    if lock=="b";
    elseif lock=="c" || lock=="m"; y = y + h*(1-factor)/2;
    elseif lock=="t"; y = y + h*(1-factor);
    else error("I don't know lock type ", lock)
    end

    h = h*factor;
    ax.set_position([x, y, w, h])

    return ax
end


"""
   ax = axisMove(xd, yd; ax=nothing)

   Moves the current axis by the indicated x and y amounts, in
   normalized figure units.

   If ax is passed then the indicated axis is modified, not necessarily
   the current one.
"""
function axisMove(xd, yd; ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax.get_position().bounds

    x += xd
    y += yd

    ax.set_position([x, y, w, h])
    return ax
end



"""
    plot_quiver(xdot, ydot; ngrid=20)

    Given functions xdot and ydot, will plot ngrid-by-ngrid quiver arrows
    within the axlim limits.
"""
function plot_quiver(xdot, ydot; ngrid=20)
    basis = (axlim[2]-axlim[1])*(0:(ngrid-1))/(ngrid-1) .+ axlim[1]
    X = repeat(basis', length(basis), 1)
    Y = repeat(basis, 1, length(basis))
    X = X[:]; Y = Y[:]

    U = xdot(X,Y)
    V = ydot(X,Y)

    quiver(X, Y, U, V, color="gray")
end

###########################################
#
#   CLINE_FOLLOWER and FIND_FIXED_POINT
#
###########################################


"""
cline_follower(func, x0, y0, xlims, ylims; Nmax=2000,
        deriv_delta=0.0001, step_delta=0.01, dxsign=1, dysign=1,
        last_branch="x")

Automatically find isoclines of two-dimensional functions.

# PARAMETERS:

- func   function that takes in x,y returns a scalar.

- x0, y0  Initial point of isocline. Whatever value that
         feval('fname', x0, y0, fpars) gives will be the z-value
         for which the isocline is computed.

- xlims   2-d double vector, specifying the smallest value of x
         and the largest value of x in which to explore the
         isocline. E.g., [0 10]. If the isocline goes outside of
         these values, the exploration is stopped.

- ylims   Bounding box as xlims (see above), but for y.


# RETURNS:

- [x, y]  Vectors such that plot(x, y) plots the isocline.


# OPTIONAL PARAMETERS:

- Nmax    Maximum number of points to explore the isocline for.

- deriv_delta   The derivatives of fname are found numerically
         by computing the function at two nearby
         points. deriv_delta defines the fraction of the bounding
         box used to identify a "nearby" point. Thus the
         x-derivative of the function is found by computing the
         function at two points separated horizontally by
         deriv_delta*diff(xlims); similarly for y.


- step_delta    Spacing between points in numerically computed
         isocline, expressed as a fraction of the bounding box.

- last_branch   Must be either the string 'x' or the string
         'y'. This is the initial direction in which the isocline
         will be explored.

- dxsign  Initial sign of horizontal direction in which to explore
         the isocline. If 'last_branch' is passed as 'y', this
         parameter is irrelevant.

- dysign  Initial sign of vertical direction in which to explore
         the isocline. If 'last_branch' is passed as 'x', this
         parameter is irrelevant.





"""
function cline_follower(func, x0, y0, xlims, ylims; Nmax=4000,
        deriv_delta=0.0001, step_delta=0.01, dxsign=1, dysign=1,
        last_branch="x")

    dx = deriv_delta*diff(xlims)[1];
    dy = deriv_delta*diff(ylims)[1];
    dydx = -1; dxdy = -1;   # Arbitrary starting points.

    step_dx = step_delta*diff(xlims)[1];
    step_dy = step_delta*diff(xlims)[1];

    x = [x0]; y = [y0]; f0 = func(x0, y0);
    if last_branch == "x"
        branch = [1];
    else
        branch = [0];
    end;

    for i=1:Nmax,
        fnow = func(x[end], y[end]);

        dfdx = (func(x[end]+dx, y[end]   ) - fnow) / dx;
        dfdy = (func(x[end]   , y[end]+dy) - fnow) / dy;
        if  abs(dfdx) <= abs(dfdy)
            # More horiz than vertical -- dydx very numerically stable
            if last_branch != "x"
                dxsign = sign(dysign*dxdy);
                last_branch = "x";
            end;
            dydx = - dfdx / dfdy;
            newx = x[end] + dxsign*step_dx;
            newy = y[end] + dxsign*step_dx*dydx;

            if (newx > xlims[2]) | (newx < xlims[1]) |
                (newy > ylims[2]) | (newy < ylims[1])
                break;
            end;

            # Now an error-correcting adjustment
            fnow = func(newx, newy);
            dfdy = (func(newx, newy+dy) - fnow) / dy;
            newy = newy + (f0 - fnow) / dfdy;
            branch = [branch ; 1];
        else
            # More vertical -- let's use dx/dy instead of dy/dx
            if last_branch != "y"
                dysign = sign(dxsign*dydx);
                last_branch = "y";
            end;
            dxdy = - dfdy / dfdx;
            newy = y[end] + dysign*step_dy;
            newx = x[end] + dysign*step_dy*dxdy;

            if (newx > xlims[2]) | (newx < xlims[1]) |
                (newy > ylims[2]) | (newy < ylims[1])
                break;
            end;

            # Now an error-correcting adjustment
            fnow = func(newx, newy);
            dfdx = (func(newx+dx, newy) - fnow) / dx;
            newx = newx + (f0 - fnow) / dfdx;
            branch = [branch ; 0];
        end;

        x = [x ; newx]; y = [y ; newy];
        if rem(i, 1000)==0; println(i); end
    end;

    return x,y,branch

end



"""
set1_index, set2_index = smallest_diff(set1, set2; midline_tol=0.3)

Given two sets of (x,y) points, finds the point in set 1 and the point in
set 2 that have the smallest Euclidean distance between them. Each set of points
should have the x coords in the first column, the y coords in the second
(i.e., should be N-by-2).

All points considered should be withing midline_tol of the x=y midline.

"""
function smallest_diff(xnull, ynull; midline_tol=0.3)
    xdiff = repeat(xnull[:,1], 1, size(ynull,1)) -
        repeat(ynull[:,1]', size(xnull,1));
    ydiff = repeat(xnull[:,2], 1, size(ynull,1)) -
        repeat(ynull[:,2]', size(xnull,1));

    alldiffs = xdiff.^2 + ydiff.^2
    sol_x = sol_y = 0
    mindiff = Inf
    for ypt=1:size(alldiffs,2)
        for xpt=1:size(alldiffs,1)
            if abs(xnull[xpt,1]-xnull[xpt,2])<midline_tol &&
                abs(ynull[ypt,1]-ynull[ypt,2])<midline_tol &&
                alldiffs[xpt,ypt] < mindiff
                mindiff = alldiffs[xpt,ypt]
                sol_x = xpt; sol_y = ypt;
            end
        end
    end

    return sol_x, sol_y
end


"""
r1, r2 = find_fp(xdot, ydot, start_pts=[0.001 0.999 ; 0.99 0.01], tol=1e-4;
    dt=0.1)

Returns two fixed points for the xdot ydot dynamics. xdot and ydot must be
functions that take x and y as variables.

start_pts has two starting points, one on each row. The function follows the
dynamics from each of those, presumably ending at two fixed points. If a change
in position after a timestep has all components smaller than tol, then the
point is considered to have reached a fixed point.

Dynamics are assumed to be confined to the positive plane: any motion component
into negative quadrants is zeroed out.
"""
function find_fp(xdot, ydot, start_pts=[0.001 0.999 ; 0.99 0.01], tol=1e-4;
    dt=0.1)

    r1  = start_pts[1,:]';
    r2  = start_pts[2,:]';

    reached_fp = false;

    while !reached_fp
        dr1dt = [xdot(r1[end,1], r1[end,2]) ydot(r1[end,1], r1[end,2])];
        newr1 = r1[end,:]' + dt*dr1dt;
        newr1[newr1.<0] .= 0;

        r1 = [r1 ; newr1];

        dr2dt = [xdot(r2[end,1], r2[end,2]) ydot(r2[end,1], r2[end,2])];
        newr2 = r2[end,:]' + dt*dr2dt;
        newr2[newr2.<0] .= 0;

        r2 = [r2 ; newr2];

        if max(maximum(abs.(r1[end,:] - r1[end-1,:])), maximum(abs.(r2[end,:] - r2[end-1,:]))) < tol
          reached_fp = true;
       end;
    end;

    return r1, r2
end


###########################################
#
#   MAIN DYNAMICS AND GRAPHICS DEFINITIONS
#
###########################################

# KEY GLOBAL PARAMETERS:
#
#  i_x      current to the x unit
#  i_y      current to the y unit
#  w_self   Self-connection weight
#  w_other  Connection to other unit
#  sigma    noise magnitude
#  dt       timestep
#  thresh   If either x or y are bigger than thresh, a decision has been reached
#


# =======================
#
#   Now the dynamics display
#
# =======================


figure(2); clf();


"""
    M, x_index, y_index, xnull, ynull = fig2_start(;quiver_only=false)

    Starts up figure 2, the dynamics plot figure.
    Opens the main dynamics axis, the two ancillary reaction time axes,
    initializes the quiver plot, threshold lines, and the tails of the dots

    Returns Jacobian for linearized dynamics near midline, the x nullcline
    (has two columns for x and y coordinates of the nullcline curve), the
    y nullcline, and the index (i.e. row) into xnull and into ynull that
    has them crossing (i.e. a fixed point) near the midline.
"""

function fig2_start(;quiver_only=false)
    global rt1ax, rt2ax, mainax, tx

    figure(2); clf(); mainax=gca();
    axisHeightChange(0.9, lock="t"); axisWidthChange(0.7, lock="r")
    axisMove(-0.1, 0.07)
    # rt1ax and rt2ax are axes to plot histograms of reaction times
    # for the x and y decisions, respectively:
    rt1ax = PyPlot.axes([0.2, 0.09, 0.55, 0.07]); axisMove(0.05, 0)
    rt2ax = PyPlot.axes([0.04, 0.25, 0.07, 0.7]); axisMove(0.05, 0)

    sca(mainax);
    # threshold lines:
    thlines = plot([0 thresh ; 1 thresh], [thresh 0 ; thresh 1],
        "-", color=[0.7, 0.7, 0.9], linewidth=1.5)
    # quiver plot handles
    q[1] = plot_quiver(xdot, ydot)
    # initialize the tails, but they are 0 length at this point
    tails[1] = plot(zeros(0,ndots), zeros(0,ndots), "-",
        color=[0.6, 0.8, 0.6], linewidth=0.5)
    if ~quiver_only
        # plot the actual dots
        p[1] = plot(xpos[:,i], ypos[:,i], "g.", markersize=12)[1]
        # if tails is true in the gui, then plot the tails
        if showtails
            # set x and y tail data for each dot:
            for k=1:ndots
                tails[1][k].set_xdata(xpos[k,1:i])
                tails[1][k].set_ydata(ypos[k,1:i])
            end
        end
    end
    xlabel("r1"); ylabel("r2")
    xlim(axlim); ylim(axlim);

    # Now for the nullclines
    r1, r2 = find_fp(xdot, ydot)  # this gives us two fixed points to explore
    # first follow the first one in two directions for the xdot nullcline:
    x1, y1 = cline_follower(xdot, r1[end,1], r1[end,2], [0,1], [0,1])
    x2, y2 = cline_follower(xdot, r1[end,1], r1[end,2], [0,1], [0,1], dxsign=-1)
    plot(x1, y1, "r"); plot(x2, y2, "r")
    xnull_x = [x1 ; x2]; xnull_y = [y1 ; y2]  # These are the x and y coords for the nullcline for x
    # then follow the other in two directions
    x1, y1 = cline_follower(xdot, r2[end,1], r2[end,2], [0,1], [0,1])
    x2, y2 = cline_follower(xdot, r2[end,1], r2[end,2], [0,1], [0,1], dxsign=-1)
    xnull_x = [xnull_x ; x1 ; x2]; xnull_y = [xnull_y ; y1 ; y2]
    plot(x1, y1, "r"); plot(x2, y2, "r")

    # now follow the first one in two directions for the ydot nullcline:
    x1, y1 = cline_follower(ydot, r1[end,1], r1[end,2], [0,1], [0,1], last_branch="y")
    x2, y2 = cline_follower(ydot, r1[end,1], r1[end,2], [0,1], [0,1], last_branch="y", dysign=-1)
    plot(x1, y1, "b"); plot(x2, y2, "b")
    ynull_x = [x1 ; x2]; ynull_y = [y1 ; y2]
    # then follow the other in two directions
    x1, y1 = cline_follower(ydot, r2[end,1], r2[end,2], [0,1], [0,1], last_branch="y")
    x2, y2 = cline_follower(ydot, r2[end,1], r2[end,2], [0,1], [0,1], last_branch="y", dysign=-1)
    ynull_x = [ynull_x ; x1 ; x2]; ynull_y = [ynull_y ; y1 ; y2]
    plot(x1, y1, "b"); plot(x2, y2, "b")


    # Collect it all so that xnull has two columns, for the x and y values
    # of the nullcline; same for ynull
    xnull = [xnull_x xnull_y]; ynull = [ynull_x ynull_y]

    # --- show evals of linearized dynamics near midline
    #
    # Let's find a crossing between the two nullclines that is near
    # the x=y midline:
    x_index, y_index = smallest_diff(xnull, ynull, midline_tol=0.3)
    x = mean([xnull[x_index,1], ynull[y_index,1]])
    y = mean([xnull[x_index,2], ynull[y_index,2]])

    # Find the Jacobian numerically for that fixed point
    eps = 1e-4
    M = [(xdot(x+eps,y)-xdot(x,y))/eps (xdot(x,y+eps)-xdot(x,y))/eps ;
         (ydot(x+eps,y)-ydot(x,y))/eps (ydot(x,y+eps)-ydot(x,y))/eps]

    # show its evals on plot
    # es = eigen(M)
    # L = es.values; V = es.vectors
    # L, V = eig(M)
    # title(@sprintf("lambda = %.3f,  %.3f", real(L[1]), real(L[2])))
    #
    # ------end of show evals

    # show current time in simulation
    #tx = text(0.8, ylim()[2]+0.05, @sprintf("t=%.3f", t),
    #    horizontalalignment="left", fontsize=12)

    return M, x_index, y_index, xnull, ynull
end


# KEY GLOBAL PARAMETERS:
#
#  i_x      current to the x unit
#  i_y      current to the y unit
#  w_self   Self-connection weight
#  w_other  Connection to other unit
#  sigma    noise magnitude
#  dt       timestep
#  thresh   If either x or y are bigger than thresh, a decision has been reached
#

"""
xdot(x, y; sigma=0)

Returns xdot given x, y, and optionally sigma. Can take x and y vectors
in which case they have to be of the same size, and the resulting xdot
is also of the same size.
"""
function xdot(x, y; sigma=0)
    dmedt = -x .+ (0.5 .+ 0.5*tanh.(i_x .+ w_self*x + w_other*y))
    if length(x)>1
        dmedt += sigma*randn(size(x))*sqrt(dt)
    else
        dmedt += sigma*randn()*sqrt(dt)
    end
    u = (x .<= 0) .& (dmedt .< 0)
    if typeof(u)<:BitArray
        dmedt[u] .= 0
    elseif u
        dmedt = 0
    end
    return dmedt
end


"""
    ydot(x, y; sigma=0)

    Returns ydot given x, y, and optionally sigma. Can take x and y vectors
    in which case they have to be of the same size, and the resulting xdot
    is also of the same size.
"""
function ydot(x, y; sigma=0)
    dmedt = -y .+ (0.5 .+ 0.5*tanh.(i_y .+ w_self*y .+ w_other*x))
    if length(y)>1
        dmedt += sigma*randn(size(y))*sqrt(dt)
    else
        dmedt += sigma*randn()*sqrt(dt)
    end
    u = (y .<= 0) .& (dmedt .< 0)
    if typeof(u)<:BitArray
        dmedt[u] .= 0
    elseif u
        dmedt = 0
    end
    return dmedt
end

###########################################
#
#   RUN MAIN DYNAMICS
#
###########################################


"""
decisionMake(;mu0=2, c=0, wR=0, wI=6.8, noise_sigma=0.05, deltat=0.05,
    bound=0.97, maxiter=100, numdots=80, displaytails=true)

Runs decision-making dynamics in a bistable system, displays plots
to illustrate trajectories and reaction times, and returns average decisions
and reaction times.

# OPTIONAL PARAMETERS:

- mu0 constant excitation to both units

- c   random dots coherence

- wR  strength of self-excitation connection

- wI  strength of inhibition connection between the two unit

- noise_sigma  standard deviation of noise added at each timestep

- deltat   timestep

- bound     value of r1 and r2 at which a decisions is deemed to be emitted,
            because of being very close to one of the stable points

- maxiter   Maximum number of timesteps in each simulation

- numdots   Number of trials to simulate in parallel

- displaytails  if true, show a green line tracing each trajectory, if false show
            only endpoints.

# RETURNS:

- p1  the fraction of decisions where unit 1 is the winner

- nd  the fraction of simulations that don't reach a decision bound.

- rt1  the mean reaction time of unit 1-winning decisions

- rt2  the mean reaction time of unit 2-winning decisions

"""
function decisionMake(;mu0=2, c=0, wR=0, wI=6.8, noise_sigma=0.1, deltat=0.05,
    bound=0.97, maxiter=150, numdots=100, displaytails=true)

    global Nmax, ndots, axlim, i, t, xpos, ypos
    global q, p, tails, thlines, mainax, rt1ax, rt2ax, tx, showtails
    global E, i_x, i_y, w_self, w_other, sigma, dt, rev, thresh

    Nmax = maxiter
    ndots = numdots
    showtails = displaytails

    # Key parameters. Current values lead to a single fixed point near  (0,0)
    E = mu0
    i_x = mu0 + 0.01*c;
    i_y = mu0 - 0.01*c;
    w_self = wR;
    w_other = - wI;
    sigma = noise_sigma;
    dt = deltat;
    rev=10000
    thresh = bound


    axlim = [-0.01, 1.01]   # these will be the limits for displaying the dynamics

    i=1; t=0      # i is the iteration number, t is the corresponding time

    xpos = zeros(ndots, Nmax)  # rows are points, columns are timestep
    ypos = zeros(ndots, Nmax)
    # Start randomly, uniformly distributed within axlim
    # xpos[:,1] = rand(ndots, 1)*0.05*(axlim[2]-axlim[1]) .+ axlim[1]
    # ypos[:,1] = rand(ndots, 1)*0.05*(axlim[2]-axlim[1]) .+ axlim[1]

    q       = Array{Any}(undef);  # q[1] will be handles to the quiver plot
    p       = Array{Any}(undef);  # p[1] will be handles to the points
    tails   = Array{Any}(undef);  # tails[1] will be handles to the tails
    thlines = Array{Any}(undef)   #
    mainax = rt1ax = rt2ax = tx = 0


    t=0; i=1
    fig2_start()
    x_rts = NaN*ones(ndots,1)
    y_rts = NaN*ones(ndots,1)

    for i=2:Nmax
        # Euler
        global rt1ax, rt2ax, all_rts
        xpos[:,i] = xpos[:,i-1] + dt*xdot(xpos[:,i-1], ypos[:,i-1]) + sigma*sqrt(dt)*randn(size(xpos,1),1)
        ypos[:,i] = ypos[:,i-1] + dt*ydot(xpos[:,i-1], ypos[:,i-1]) + sigma*sqrt(dt)*randn(size(ypos,1),1)
        t = t+dt
        # tx.set_text(@sprintf("t=%.3f", t))
        x_rts[isnan.(x_rts) .& (xpos[:,i] .> thresh) .& (xpos[:,i-1] .<= thresh)] .= t
        y_rts[isnan.(y_rts) .& (ypos[:,i] .> thresh) .& (ypos[:,i-1] .<= thresh)] .= t
        all_rts = [x_rts[.!isnan.(x_rts)] ; y_rts[.!isnan.(y_rts)]]
    end

    if any(.!isnan.(x_rts))
        sca(rt1ax); cla();
        plt.hist(x_rts[.!isnan.(x_rts)], range(0, stop=maximum(all_rts), length=20))
        xlabel("reaction time")
    end
    if any(.!isnan.(y_rts))
        sca(rt2ax); cla();
        plt.hist(y_rts[.!isnan.(y_rts)],
            range(0, stop=maximum(all_rts), length=20), orientation=:horizontal)
        ylabel("reaction time")
    end

    sca(mainax)

    p[1].set_xdata(xpos[:,Nmax])
    p[1].set_ydata(ypos[:,Nmax])
    if showtails
        for k=1:ndots
            tails[1][k].set_xdata(xpos[k,1:Nmax])
            tails[1][k].set_ydata(ypos[k,1:Nmax])
        end
    end

    r1s = x_rts[findall(.!isnan.(x_rts))]
    r2s = y_rts[findall(.!isnan.(y_rts))]
    ndecs = length(r1s)+length(r2s)

    println(@sprintf("P1=%.2f, mean RT1=%.2f, mean RT2=%.2f",
        length(r1s)/ndecs, mean(r1s), mean(r2s)))
    return length(r1s)/ndecs, (numdots - ndecs)/ndecs, mean(r1s), mean(r2s)
end

end  # end Module BistableDynamics

using .BistableDynamics
using PyPlot
pygui(true)
