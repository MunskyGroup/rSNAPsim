# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 14:17:56 2020

@author: willi
"""


    def elongation_animation(self, ti=0, tf=1000, tstep=1000, cell_radius=50, imagesize=5, dpi=90, filename='simulated_cell', ssa_obj=None, fcolor='#00FF00' ,rnacolor='#FF0000', xkcd=False):
        '''
        function that creates a mrna translation animation
        '''

        custom_cmap = ['#69dd42', '#e5361b', '#db11c7']
        def rpts(x, y, angle):
            nx = np.cos(angle)*x - np.sin(angle)*y
            ny = np.sin(angle)*x + np.cos(angle)*y
            return nx, ny


        def update_line(num, xpos, ypos, line):  #function for the FuncAnimation
            if num != 0:
                ax.get_lines()[-1].remove()

                for child in ax.get_children():  #remove the previous patch collection (green spots)

                    if isinstance(child, PatchCollection):
                        child.remove()


            patches = []
            gp = []
            ep = []
            radi = np.ones(xpos[:, inds[num]].shape)*4 #create a max radius of 3 for intensity vecs
            ypos = np.ones(xpos[:, inds[num]].shape)*(ytop+3)
            x = xpos[:, inds[num]]
            x[np.where(x == 0)] = x[np.where(x == 0)] - 300

            for x1, y1, r in zip(xpos[:, inds[num]], ypos, radi):   #make circle objects of radius based on ivec
                circle = mpatches.Circle((x1, y1), r, facecolor='#FF0000', edgecolor='k')
                patches.append(circle)

            pcolor = custom_cmap[0]

            for i in range(len(x.flatten())):

                if x[i] > 0:
                    xpts = np.linspace(0, int(x[i])-1, int(x[i]))
                    ypts = 5*np.sin(1/10*np.linspace(0, int(x[i])-1, int(x[i])))
                    xpts, ypts = rpts(ypts, xpts, 1)
                    ypts = ypts+ytop+3
                    xpts = xpts+x[i]
                    radi = np.ones(xpts.shape)*1
                    k = 0
                    ypts = np.fliplr(np.atleast_2d(ypts))
                    ypts = ypts.flatten()
                    xpts = np.fliplr(np.atleast_2d(xpts))
                    xpts = xpts.flatten()

                    for x2, y2, r2 in zip(xpts, ypts, radi):
                        probloc = False
                        j = 0
                        for key in epitopes.keys():
                            if k in epitopes[key]:
                                probloc = True
                                pcolor = custom_cmap[j]
                                j += 1
                        rx = np.random.rand()*2
                        ry = np.random.rand()*2
                        if probloc == False:

                            circle = mpatches.Circle((x2+rx, y2+ry), r2, facecolor='#0000FF', edgecolor='#FFFFFF', lw=2, ls='solid')
                            gp.append(circle)
                        else:
                            circle = mpatches.Circle((x2+rx, y2+ry), r2*3, facecolor='#00FF00', edgecolor='#000000', lw=2, ls='solid')
                            ep.append(circle)

                        k += 1


                #fig.gca().add_artist(circle)
            '''
            xs = np.flip(np.sort(xpos[:,inds[num]][0].flatten()),axis=0)
            for i in range(max_ribs):
                line.set_data(xpos[:,inds[num]],ypos[inds[num]])
                line.set_linewidth(0)
                line.set_marker('o')
                line.set_markersize(3)
            '''
            p = PatchCollection(patches, facecolors=('#FF0000',), zorder=5)  #create a patch collection to add to axis

            m = PatchCollection(gp, facecolors=('#0000FF',), lw=2, zorder=3)  #create a patch collection to add to axis
            e = PatchCollection(ep, facecolors=(pcolor,), zorder=4)



            n = num

            ax.plot(np.linspace(0, tag_length, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]))[:n], 3*ssa_obj.intensity_vec.flatten()[:n]+total_length, color=pcolor)

            fldot = mpatches.Ellipse((total_length-30, total_length+40), width=ssa_obj.intensity_vec.flatten()[n], height=ssa_obj.intensity_vec.flatten()[n]*1.0, color=pcolor)
            f = [fldot]
            fe = PatchCollection(f, facecolors=(pcolor,), zorder=4)
            ax.add_collection(p)  #adds the circles to axis
            ax.add_collection(m)  #adds the circles to axis
            ax.add_collection(e)
            ax.add_collection(fe)
            plt.xlabel(str(inds[num]))  #update time label
            return line,

        if ssa_obj == None:
            ssa_obj = self.ssa_solver(n_traj=1, tf=tf, tstep=tstep)
        if xkcd == True:
            plt.xkcd()



        fig1 = plt.figure(figsize=(imagesize+5, imagesize), dpi=dpi)  #make figure
        fig1.tight_layout()

        ax = fig1.add_subplot('111')
        ax.set_aspect(1)

        tag_length = self.POI.tag_length
        total_length = self.POI.total_length

        epitopes = self.POI.tag_epitopes

        tag_length = total_length - self.POI.gene_length


        ax.cla()
        ybot = 90
        ytop = 110
        ax.plot([0, total_length], [ybot, ybot], color='white', zorder=3)
        ax.plot([0, total_length], [ytop, ytop], color='white', zorder=3)
        ax.plot([0, 0], [ybot, ytop], color='white', zorder=3)
        ax.plot([total_length, total_length], [ybot, ytop], color='white', zorder=3)
        ax.axis([-10, total_length+10, 80, total_length+np.max(ssa_obj.intensity_vec)*3+20])
        ax.plot([tag_length, tag_length], [ybot, ytop], color='white', linewidth=1, zorder=3)
        k = 0
        for key in epitopes.keys():
            for i in range(len(epitopes[key])):
                ax.plot([epitopes[key][i], epitopes[key][i]], [ybot, ytop], color=custom_cmap[k], linewidth=2, zorder=3)
            rect = mpatches.Rectangle(xy=(tag_length, ybot), width=total_length-tag_length, height=ytop-ybot, color='#0000FF')
            #ax.fill_between([tag_length,tag_length,total_length,total_length],[ybot,ytop,ytop,ybot],color='#00FF00')
            ax.add_patch(rect)
            k += 1

        ticks = np.linspace(0, total_length, 10).astype(int)
        ax.set_xticks(ticks)
        ax.set_xlabel('Codon Position')
        ax.get_yaxis().set_visible(False)
        ax.set_facecolor('k')
        filename = 'elong.gif'

        Writer = animation.writers['pillow']

        print('making movie...')
        max_ribs = np.max(np.nonzero(ssa_obj.solutions[0])[0])

        l, = plt.plot([], [], 'r-')
        t = ssa_obj.time_vec_fixed[ssa_obj.start_time:]
        inds = np.linspace(0, len(t)-1, len(t)).astype(int)
        xpos = np.zeros((max_ribs, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:])))

        ypos = np.ones((1, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]))).flatten()

        xpos[:, :] = ssa_obj.solutions[0][:max_ribs, ssa_obj.start_time:len(ssa_obj.time_vec_fixed)]

        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
        line_ani = animation.FuncAnimation(fig1, update_line, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]), fargs=(xpos, ypos, l),
                                           interval=50, blit=True)

        line_ani.save((filename), writer=writer)  #save the animation




    def simulate_cell(self, diffusion_constant, kon, koff, kRNA, kdecay, ti=0, tf=1000, tstep=1000, cell_radius=50, imagesize=5, dpi=90, filename='simulated_cell', ssa_obj=None, fcolor='#00FF00', rnacolor='#FF0000'):
        '''
        [DNA] ==kRNA==> [RNA] <==koff== [RNA*] ==translation simulation==> [Protein]===> null
                        // ||             /\
                        || `'=====kon====='`
                        ||
                        \/
                        null

        '''
        print('simulating RNA creation....')
        t = np.linspace(ti, tf, tstep)


        dna_s = np.array([[ 0,  0],
                          [ 1, -1]])

        dna_w1 = np.array([[kRNA, 0],
                           [0, 0]],dtype=float)


        dna_w0 = np.array([[0], [0]])


        dna_si = GenericSSA(type='linear' )
        dna_si.W1 = dna_w1
        dna_si.W0 = dna_w0
        dna_si.S = dna_s

        dna_si.ti = t[0]
        dna_si.tf = t[-1]
        dna_si.n = 1
        xi = np.zeros((2, 1))
        xi[0] = 1
        dna_si.xi = xi
        dna_si.ptimes = len(t)

        dna_si.time_variant = False
        dna_si._solve(1)
        rna_creation_data = dna_si.data




        stoich = np.array([[  0,    0,  1],
                           [  -1,  1, -1],
                           [  1, -1, 0]])

        propensity = np.array([[0, kon, 0],
                              [0, 0,koff],
                              [0,kdecay, 0]], dtype=float)

        w0 = np.array([[0],[0],[0]])

        solver_instance = GenericSSA(type='linear' )
        solver_instance.W1 = propensity
        solver_instance.W0 = w0
        solver_instance.S = stoich

        solver_instance.ti = t[0]
        solver_instance.tf = t[-1]
        solver_instance.n = 1
        xi = np.zeros((3,1))
        xi[1] = 1
        solver_instance.xi = xi
        solver_instance.ptimes = len(t)

        solver_instance.time_variant = False




        print('simulating RNA activation....')



        R = cell_radius
        squarelen = float(R/np.sqrt(2))

        n_RNA_t = np.zeros((len(t),int(np.max(rna_creation_data[1]))))

        nRNA = 0
        nparticles = (int(np.max(rna_creation_data[1])))
        for i in range(len(t)):

            while nRNA != rna_creation_data[1][i]:
                data = solver_instance._solve(1)


                rnaonoff = data[2] + 1 - data[0]



                n_RNA_t[i:, nRNA] = rnaonoff[:-i].flatten()
                nRNA += 1


        rna_particles = n_RNA_t.T
        rna_exist = np.where(rna_particles >0,1,0)
        rnaex = data



        print('simulating RNA motion....')
        rna_locations = np.empty((nparticles, len(t), 2))

        dt = t[-1]/len(t)

        delta = diffusion_constant



        def linecirc(m, b, xc, yc, r):

            if np.isinf(m) == False:
                a = 1+m**2
                e = 2*(m*(b-yc)-xc)
                c = yc**2+xc**2 + b**2-2*yc*b-r**2
                x = np.roots([a, e, c])

                if np.isreal(x).all() == False:
                    x = [np.nan, np.nan]
                    y = [np.nan, np.nan]
                else:
                    y = [b + m*x[0], b+m*x[1]]

            elif abs(xc-b) > r:
                x = [np.nan, np.nan]
            else:
                x = [b, b]
                step = np.sqrt(r**2-(b-xc)**2)
                y = [yc + step, yc-step]

            return [x[0], y[0]], [x[1], y[1]]


        def dist(x1, y1, x2, y2):
            return np.sqrt((x1-x2)**2+(y1-y2)**2)


        for i in range(nparticles):
            x = np.empty((2,len(t) - np.where(rna_exist[i] != 0 )[0][0]  ))
            centers = np.zeros(x.shape)
            x[:,0] = np.random.random()*squarelen
            x0 = [  ((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4),((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4) ]


            x0 = x0 - np.array([R, R])
            x[:,0] =x0
            r = norm.rvs(size=np.array(x0).shape + (len(t) - np.where(rna_exist[i] !=0 )[0][0],), scale=delta*np.sqrt(dt))



            out = np.empty(r.shape)

            np.cumsum(r, axis=-1, out=out)
            out += np.expand_dims(np.array(x0), axis=-1)

            #out = np.array([[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36],
                            #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]])

            centers = np.zeros(out.shape)
            dists = np.zeros((x.shape[1], 1)).flatten()

            incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
            dists[np.where(out[0] != 0)] = incirc[np.where(out[0] != 0)]

            while len(np.where(dists>R)[0]) != 0:   #trajectory left the cell
                out = out.T
                left_cell = np.where(dists > R)[0][0]



                pts = [[out[left_cell][0], out[left_cell][1]], [out[left_cell-1][0], out[left_cell-1][1]]]

                p = np.polyfit([out[left_cell][0], out[left_cell-1][0]], [out[left_cell][1], out[left_cell-1][1]], 1)
                m = p[0]
                b = p[1]

                intercepts = linecirc(m, b, 0, 0, R)
                if dist(*tuple(intercepts[0])+tuple(pts[0])) > dist(*tuple(intercepts[1])+tuple(pts[0])):
                    inter = np.array(intercepts[1])
                else:
                    inter = np.array(intercepts[0])

                a = out[left_cell] - inter


                out[left_cell-1:] = out[left_cell-1:] - 2*(np.dot(inter, a)/np.linalg.norm(inter)**2)*inter




                dists = np.zeros((x.shape[1], 1)).flatten()
                out = out.T
                incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
                dists[np.where(out[0] != 0)] = incirc[np.where(out[0] != 0)]



            data = ((out.T).T*rna_exist[i][np.where(rna_exist[i] != 0)[0][0]:].T).T
            data[np.where(rna_exist[i] != 0)[0][-1]- np.where(rna_exist[i] != 0)[0][0]+1 :] = -R


            rna_locations[i, np.where(rna_exist[i] != 0)[0][0]:, :] =  data
            rna_locations[i, :np.where(rna_exist[i] != 0)[0][0], :] =  -R


        print(nparticles)
        rna_loc_compressed = rna_locations[np.where(np.sum(np.sum(rna_locations+R, axis=1), axis=1) > 0)]

        if ssa_obj == None:
            print('no ssa data given')
            print('simulating translation....')
            print(int(rna_loc_compressed.shape[0]))
            ssa_obj = self.ssa_solver(n_traj=int(rna_loc_compressed.shape[0]),tf=tf,tstep=tstep)

            ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
            ivec = ivec.T  #get the intensity vec for the "fluorescence"


        else:
            print('Translation data given')
            print('Given ' + str(ssa_obj.n_traj) + ' Needed '+str(int(rna_loc_compressed.shape[0])) )
            if ssa_obj.n_traj  < int(rna_loc_compressed.shape[0]):
                print('simulating ' + str(int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj) + ' additional trajectories....')
                ssa_obj = self.ssa_solver_append(ssa_obj, n=int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj)
                ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
                ivec = ivec.T  #get the intensity vec for the "fluorescence"

            else:
                ivec = ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])]/np.max(ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])])
                ivec = ivec[0:int(rna_loc_compressed.shape[0])].T  #get the intensity vec for the "fluorescence"







        print('making movie...')
        #simulate brownian motion
        def update_line(num, xpos,ypos, line):  #function for the FuncAnimation
            if num !=0:

                for child in ax.get_children():  #remove the previous patch collection (green spots)

                    if isinstance(child, PatchCollection):
                        child.remove()
                    if isinstance(child, mpatches.Ellipse):
                        child.remove()

            patches = []
            radi = 3*ivec[inds[num]]   #create a max radius of 3 for intensity vecs


            for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
                circle = mpatches.Circle((x1, y1), r, color=fcolor)
                patches.append(circle)
                #fig.gca().add_artist(circle)

            line.set_data(xpos[inds[num]], ypos[inds[num]])
            line.set_linewidth(0)
            line.set_marker('o')
            line.set_markersize(3)
            line.set_color(rnacolor)
            line.set
            p = PatchCollection(patches, zorder=3, facecolors=(fcolor,))  #create a patch collection to add to axis
            ax.add_collection(p)  #adds the circles to axis


            p = mpatches.Circle((0,0), radius=R, color='black')  #add the black circle
            ax.add_patch(p)


            whitep = mpatches.Ellipse((-R, -R), width=7, height=7, color='white', zorder=5)  #add the black circle
            ax.add_patch(whitep)

            plt.xlabel(str(inds[num]))  #update time label


            return line,
        xpos = rna_loc_compressed.T[0]
        ypos = rna_loc_compressed.T[1]


        filetype='.mov'

        if filetype == '.gif':

            Writer = animation.writers['pillow']
        if filetype == '.html':
            Writer = animation.writers['html']
        if filetype == '.gif':

            Writer = animation.writers['FFMpeg']

        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

        fig1 = plt.figure(figsize=(imagesize, imagesize),dpi=dpi)  #make figure
        fig1.tight_layout()

        ax= fig1.add_subplot('111')
        plt.yticks([])
        plt.xticks([])
        p = mpatches.Circle((0, 0), radius=R, color='black')  #add the black circle
        ax.add_patch(p)
        plt.gca().set_aspect('equal', adjustable='box')

        l, = plt.plot([], [], 'r-')
        plt.xlim(-R-10, R+10)
        plt.ylim(-R-10, R+10)
        plt.xlabel('0')
        plt.title('Simulated Cell')

        inds = np.linspace(0, len(t)-1, len(t)).astype(int)
        #creates the animation
        line_ani = animation.FuncAnimation(fig1, update_line, tstep, fargs=(xpos,ypos, l),
                                           interval=50, blit=True)

        line_ani.save((filename + filetype), writer=writer)  #save the animation



        #return solver_instance,n_RNA_t,rna_creation_data,data,rna_locations
        return rna_locations, rna_loc_compressed, rna_particles, rna_creation_data, rna_exist, rnaonoff, rnaex


    def get_simulated_mov(self, ssa_obj, filename, filetype):
        '''
        Create a gif or html file of the simulated circlular cell from any sms ssa object

        '''

        R = 50   #set up the random circle and points
        num = ssa_obj.n_traj
        r1 = np.zeros((1, num)).flatten()
        theta1 = np.zeros((1, num)).flatten()
        x1 = np.zeros((1, num)).flatten()
        y1 = np.zeros((1, num)).flatten()


        r2 = np.zeros((1, num)).flatten()
        theta2 = np.zeros((1, num)).flatten()
        x2 = np.zeros((1, num)).flatten()
        y2 = np.zeros((1, num)).flatten()

        for n in range(0, num):  #for  all trajectories make initial points
            r1[n] = R*np.sqrt(np.random.random(1))
            r2[n] = R*np.sqrt(np.random.random(1))
            theta1[n] = 2*np.pi*np.random.random(1)
            theta2[n] = 2*np.pi*np.random.random(1)
            x1[n] = np.cos(theta1[n])*r1[n]
            x2[n] = np.cos(theta2[n])*r2[n]
            y1[n] = np.sin(theta1[n])*r1[n]
            y2[n] = np.sin(theta1[n])*r2[n]

        movement = .7
        xpos = np.zeros((len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]), num))
        ypos = np.zeros((len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]), num))

        #over time perterub the simulated ribosome posistions and save to xpos , ypos
        for j in range(len(ssa_obj.time_vec_fixed[ssa_obj.start_time:])):
            if j == 0:
                xpos[0] = x1
                ypos[0] = y1
            else:
                for i in range(0,num):
                    xpos[j, i] = xpos[j-1,i]-movement + 2*movement*np.random.random(1)
                    if xpos[j, i] > 52:
                        xpos[j, i] = 51
                    if xpos[j, i] < -52:
                        xpos[j, i] = -51
                    ypos[j, i] = ypos[j-1,i]-movement + 2*movement*np.random.random(1)
                    if ypos[j, i] > 52:
                        ypos[j, i] = 51
                    if ypos[j, i] < -52:
                        ypos[j, i] = -51

        ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
        ivec = ivec.T  #get the intensity vec for the "fluorescence"
        k = 0
        def update_line(num, xpos, ypos, line):  #function for the FuncAnimation
            if num !=0:
                for child in ax.get_children():  #remove the previous patch collection (green spots)

                    if isinstance(child, PatchCollection):
                        child.remove()

            patches = []
            radi = 3*ivec[inds[num]]   #create a max radius of 3 for intensity vecs


            for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
                circle = mpatches.Circle((x1, y1), r, color='#00FF00')
                patches.append(circle)
                #fig.gca().add_artist(circle)

            line.set_data(xpos[inds[num]],ypos[inds[num]])
            line.set_linewidth(0)
            line.set_marker('o')
            line.set_markersize(.5)
            p = PatchCollection(patches, zorder=2, facecolors=('#00FF00',))  #create a patch collection to add to axis





            ax.add_collection(p)  #adds the circles to axis
            plt.xlabel(str(inds[num]))  #update time label





            return line,
        if filetype == '.gif':

            Writer = animation.writers['pillow']
        if filetype == '.html':
            Writer = animation.writers['html']
        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

        fig1 = plt.figure()  #make figure
        fig1.tight_layout()

        ax= fig1.add_subplot('111')
        plt.yticks([])
        plt.xticks([])
        p = mpatches.Circle((0, 0), radius=65, color='black')  #add the black circle
        ax.add_patch(p)
        plt.gca().set_aspect('equal', adjustable='box')

        l, = plt.plot([], [], 'r-')
        plt.xlim(-70, 70)
        plt.ylim(-70, 70)
        plt.xlabel('0')
        plt.title('Simulated Cell')
        inds = np.linspace(0, xpos.shape[0]-1, 120).astype(int)
        #creates the animation
        line_ani = animation.FuncAnimation(fig1, update_line, 120, fargs=(xpos,ypos, l),
                                           interval=50, blit=True)

        line_ani.save((filename + filetype), writer=writer)  #save the animation