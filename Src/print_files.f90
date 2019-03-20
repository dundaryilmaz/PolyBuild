module print_files

  use polymer_type
  use ceramic_type
  use global_variables
  use polymer_module
  use nanoparticle_module
  implicit none
  integer, allocatable , dimension(:,:) :: bondtable
  integer :: n_total_bonds 
  contains
 
  function connection_table(filename) result(iostat)
    character(len=*) :: filename
    character(len=40) :: fname
    integer ::i , j , k,q,jatom,qatom,katom,iatom
    real*8 , dimension(3) :: vec1,vec2,vec12
    integer :: n_atoms,iostat,cnt,nbonds,nangles , ntorsions,maxtype
    real*8 , dimension(:,:) , allocatable :: coords
    character(len=2),dimension(:),allocatable :: atom_symbols
    integer, dimension(:,:),allocatable :: bonds
    integer ,dimension(20) :: temp_bonds
    real*8 :: dist
    integer,dimension(:,:),allocatable ::  bondlist,tempblist,anglelist,tempalist
    integer, dimension(:,:),allocatable :: temptlist,torsanglelist
    !---------------------------------
    ! Detect bonds
    print*,'CONNECTION TABLE'
    fname=trim(filename)//'opls.lmp'
    open(9,file=fname)
  !  open(133,file=trim(filename)//'.xyz')
  !  read(133,*) n_atoms
    print*,'Natoms : ',size(printedatoms)
    print*,'Nbonds : ', n_total_bonds
    n_atoms=size(printedatoms)
  !  if ( n_amorph_chains == 0 ) return
    allocate(coords(n_atoms,3),atom_symbols(n_atoms),bonds(n_atoms,10))
 !   read(133,*)
 !   do i=1, n_atoms
 !     read(133,*) atom_symbols(i), coords(i,:)
 !   end do
    nbonds=0
    allocate( tempblist(n_atoms*100,3))
    do i=1, n_atoms
      vec1=printedatoms(i)%coord(:)
      cnt=0
      bonds(i,1)=cnt
      !print*,'i ',i
      temp_bonds=0
      

      !if( printedatoms(i)%symbol .eq. 'H ') cycle
      do j=1, n_atoms
        !print*,'j ',j
        if ( i == j ) cycle
        if (printedatoms(i)%molid > n_boxes ) cycle
        if ( printedatoms(i)%molid /= printedatoms(j)%molid ) cycle
        vec2=printedatoms(j)%coord(:)
        dist=pbc_dist(vec1,vec2)
        if (  dist < cutoff(printedatoms(i)%symbol,printedatoms(j)%symbol) ) then
          !print*, i,j,atom_symbols(i),atom_symbols(j),dist
          !if ( dist > 1.09 ) then
          cnt=cnt+1
          temp_bonds(cnt)=j
          if( j > i ) then
            nbonds=nbonds+1
            tempblist(nbonds,1)=bondtype(printedatoms(i)%symbol,printedatoms(j)%symbol)
            tempblist(nbonds,2)=i
            tempblist(nbonds,3)=j
            write(*,'(A,x,A,x,i4,x,A,x,i4,x,A,x,i4)') &
&                'iatom : ','Chain : ', printedatoms(i)%molid, 'Mer : ', printedatoms(i)%mer_index ,&
&                'Monomer Atom Index : ', printedatoms(i)%atom_mer_index                 
            write(*,'(A,x,A,x,i4,x,A,x,i4,x,A,x,i4)') &
&                'jatom : ','Chain : ', printedatoms(j)%molid, 'Mer : ', printedatoms(j)%mer_index ,&
&                'Monomer Atom Index : ', printedatoms(j)%atom_mer_index    
            print*, 'Distance : ', dist             
          end if
      !end if 

        end if
      end do
      bonds(i,1)=cnt
      bonds(i,2:cnt+1)=temp_bonds(1:cnt)
      !print*,i,bonds(i,2:cnt+1)
    end do
    allocate( bondlist(n_total_bonds,3))
    do i=1, n_total_bonds
       iatom = bondtable(i,1)
       jatom = bondtable(i,2)
       bondlist(i,1) = bondtype(printedatoms(iatom)%symbol,printedatoms(jatom)%symbol)
       bondlist(i,2) = iatom
       bondlist(i,3) = jatom
    end do

    deallocate(tempblist)
    write(*,*) 'Bonds'
    !write(*,*)
    ! bond angle
    allocate(tempalist(n_atoms*1000,4))
    nangles=0
    do i=1,n_atoms
      if( bonds(i,1)  < 2 ) cycle
      if( printedatoms(i)%symbol .eq. 'H ') cycle
      do j=1,bonds(i,1)-1
        jatom=bonds(i,j+1)
        do k=j+1,bonds(i,1)
          katom=bonds(i,k+1)
          nangles=nangles+1
          tempalist(nangles,1)=angletype(printedatoms(jatom)%symbol,printedatoms(i)%symbol,printedatoms(katom)%symbol)
          tempalist(nangles,2)=jatom
          tempalist(nangles,3)=i
          tempalist(nangles,4)=katom
          iostat=is_bonded(bonds,i,jatom)
          iostat=is_bonded(bonds,i,katom)
          !print*,i,bonds(i,j+1),bonds(i,k+1)
        end do
      end do
    end do
    allocate(anglelist(nangles,4))
    anglelist=tempalist(1:nangles,:)
    allocate(temptlist(n_atoms*2000,5))
    !torsion angle
    ntorsions=0
    do i=1,n_atoms
      if ( printedatoms(i)%symbol .eq. 'H ') cycle
      if ( bonds(i,1) <2   ) cycle
      do j=1,bonds(i,1)
        jatom=bonds(i,1+j)
        if ( jatom < i ) cycle
        if ( printedatoms(jatom)%symbol .eq. 'H ') cycle
        if( bonds(jatom ,1) < 2 ) cycle
        do k=1,bonds(jatom,1)
          katom=bonds(jatom,1+k)
          if( katom .eq. i ) cycle
          do q=1,bonds(i,1)
            if( q .eq. j ) cycle
            qatom=bonds(i,1+q)
            ntorsions=ntorsions+1
            temptlist(ntorsions,1)=torsiontype(printedatoms(qatom)%symbol,printedatoms(i)%symbol, &
&                           printedatoms(jatom)%symbol , printedatoms(katom)%symbol)
             temptlist(ntorsions,2)=qatom
             temptlist(ntorsions,3)=i
             temptlist(ntorsions,4)=jatom
             temptlist(ntorsions,5)=katom

!            print*, qatom,i,jatom,katom,torsiontype(printedatoms(qatom)%symbol,printedatoms(i)%symbol, &
!&                           printedatoms(jatom)%symbol , printedatoms(katom)%symbol)

          end do
        end do
      end do
    end do
    print*,'Connection Table Done.'
    allocate(torsanglelist(ntorsions,5))
    torsanglelist=temptlist(1:ntorsions,:)

    maxtype=maxval(printedatoms(:)%atomtype)
!maxval(base_mer%atype(:))
    write(9,*) 
    write(9,*) n_atoms, ' atoms'
    write(9,*) maxtype , ' atom types' 
    write(9,*) n_total_bonds , ' bonds'
    write(9,*) nangles, ' angles'
    write(9,*) ntorsions , ' dihedrals'




    write(9,*) maxval(bondlist(:,1)), '  bond types'
    write(9,*) maxval(anglelist(:,1)),' angle types'
    write(9,*) maxval(torsanglelist(:,1)),' dihedral types'
    write(9,*) 0, ' improper types '
    write(9,*) 


    write(9,1005) 0.d0 , BOX(1)
    write(9,1006) 0.d0 , BOX(2)
    write(9,1007) 0.d0 , BOX(3)

    write(9,*) 
    write(9,*) ' Atoms'
    write(9,*) 
    do i=1,n_atoms

      write(9,'(i6,x,i3,x,i3,x,4(f12.6,x))') i,printedatoms(i)%molid,printedatoms(i)%atomtype,&
&           printedatoms(i)%charge,printedatoms(i)%coord
    end do

    write(9,*) 
    write(9,*) ' Bonds'
    write(9,*) 
    do i=1,n_total_bonds
      write(9,*) i,bondlist(i,:)
    end do
    write(9,*) 
    write(9,*) ' Angles'
    write(9,*) 
    do i=1,nangles
      write(9,*) i, anglelist(i,:)
    end do
    write(9,*) 
    write(9,*) ' Dihedrals'
    write(9,*) 
    do i=1,ntorsions
      write(9,*) i,torsanglelist(i,:)
    end do


        



  
1005    format(2(f16.4,2x) , 'xlo xhi')
1006    format(2(f16.4,2x) , 'ylo yhi')
1007    format(2(f16.4,2x) , 'zlo zhi')

  end function
  function is_bonded(bonds,iatom,jatom) result (iostat)
    integer , dimension (:,:) ::bonds
    integer  :: iatom, jatom,i,j
    integer :: n_bonds
    integer :: iostat,jp
    logical :: bonded
    bonded=.false.
     !print*, size(bonds,dim=1)
     do i=1, bonds(iatom,1)
         jp=bonds(iatom,i+1)
         if( jp .eq. jatom) bonded=.true.
     end do
      if ( .not. bonded) then
          print*,iatom,jatom
            print*,'iatom bonds : ',bonds(iatom,:)
     pause
      end if
     iostat=0


  end function
  function torsiontype(symL,symCL,symCR,symR)
    character(len=2) ::symL,symCL,symCR,symR 
    integer torsiontype
      torsiontype=1
      if(symL == 'C ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='C ') torsiontype=1
      if(symL == 'H ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='C ') torsiontype=2
      if(symL == 'C ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='H ') torsiontype=2
      if(symL == 'H ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='H ') torsiontype=3
      if(symL == 'O ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='C ') torsiontype=4
      if(symL == 'C ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='O ') torsiontype=4
      if(symL == 'O ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='H ') torsiontype=5
      if(symL == 'H ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='O ') torsiontype=5
      if(symL == 'O ' .and. symCL == 'C '.and. symCR =='C '.and. symR =='O ') torsiontype=6
!     write(*,'(4(A2,x),i2)') syml,symcl,symcr,symr,torsiontype
  end function torsiontype

  function angletype(symL,symC,symR)
    character(len=2) :: symL,symC,symR
    integer  angletype
      angletype=0
      if(symL == 'C ' .and. symC == 'C '.and. symR =='C ') angletype=1
      if(symL == 'H ' .and. symC == 'C '.and. symR =='C ') angletype=2
      if(symL == 'C ' .and. symC == 'C '.and. symR =='H ') angletype=2
      if(symL == 'H ' .and. symC == 'C '.and. symR =='H ') angletype=3
      if(symL == 'O ' .and. symC == 'C '.and. symR =='C ') angletype=4
      if(symL == 'C ' .and. symC == 'C '.and. symR =='O ') angletype=4
      if(symL == 'H ' .and. symC == 'C '.and. symR =='O ') angletype=5
      if(symL == 'O ' .and. symC == 'C '.and. symR =='H ') angletype=5
      if(symL == 'O ' .and. symC == 'C '.and. symR =='O ') angletype=6
      if(symL == 'O ' .and. symC == 'O '.and. symR =='C ') angletype=7
      if(symL == 'C ' .and. symC == 'O '.and. symR =='O ') angletype=7
      if (angletype .eq. 0 ) then
      !  print*,symL,symC,symR
      !   pause
end if
      
  end function

 
  function cutoff(symA,symB) 
    character(len=2) :: symA,symB
    real*8 :: cutoff

      cutoff=1.6
      if(symA == 'C ' .and. symB == 'C ') cutoff=1.7d0
      if(symA == 'C ' .and. symB == 'H ') cutoff=1.3d0
      if(symA == 'H ' .and. symB == 'C ') cutoff=1.3d0
      if(symA == 'H ' .and. symB == 'H ') cutoff=0.1d0
      if(symA == 'C ' .and. symB == 'O ') cutoff=1.6d0
      if(symA == 'O ' .and. symB == 'C ') cutoff=1.6d0
      if(symA == 'O ' .and. symB == 'H ') cutoff=1.3d0
      if(symA == 'H ' .and. symB == 'O ') cutoff=1.3d0

  end function

  function atomtype(symA) 
    character(len=2) :: symA
    integer :: atomtype

      atomtype=1
      if(symA == 'C ' ) atomtype=1
      if(symA == 'H ' ) atomtype=2
      if(symA == 'O ' ) atomtype=3
   end function

  function bondtype(symA,symB) 
    character(len=2) :: symA,symB
    integer :: bondtype

      bondtype=1
      if(symA == 'C ' .and. symB == 'C ') bondtype=1
      if(symA == 'C ' .and. symB == 'H ') bondtype=2
      if(symA == 'H ' .and. symB == 'C ') bondtype=2
      if(symA == 'H ' .and. symB == 'H ') bondtype=3
      if(symA == 'O ' .and. symB == 'C ') bondtype=4
      if(symA == 'C ' .and. symB == 'O ') bondtype=4
      if(symA == 'O ' .and. symB == 'O ') bondtype=5
      if(symA == 'O ' .and. symB == 'H ') bondtype=6
      ! print*,symA,symB,bondtype

  end function


  function print_all( chains, blocks , filename,blocknumber,particle) result (iostat)
    type(chain)   ,dimension(:) :: chains
    type(ceramic) ,dimension(:) :: blocks
    type(nanoparticle) , optional :: particle
    integer , optional :: blocknumber
    integer :: nblocks
    character(len=*) :: filename
    character(len=50) ::filexyz, filelammps,filecfg,filebgf
    integer :: iostat,mer_id
    integer :: n_total_atoms,n_polymer_atoms,n_ceramic_atoms,n_nanopartatoms
    integer :: i , j , k,cnt,temp,maxtype,maxtypeo,q
    integer :: cbond,hc,tc,iatom, jatom, b
    real*8 ,dimension(3) :: rvec
    logical :: radical
    nblocks=size(blocks)
    if( present(blocknumber) ) nblocks=blocknumber
!   print*,'-',nblocks
    filelammps=trim(filename)//'.lmp'
    filecfg=trim(filename)//'.cfg'
    filexyz=trim(filename)//'.xyz'
    filebgf=trim(filename)//'.bgf'
 !  print*,'--'
    temp=0
    n_ceramic_atoms=0
    maxtype=1
    if( nblocks .gt.0 ) then
      n_ceramic_atoms=sum(blocks(:)%n_atoms)
      do i=1,nblocks
        temp=temp+sum(blocks(i)%deleted(:))
        if( sum(blocks(i)%deleted(:)) > 0 ) &
          write(*,'(A,i3,x,i6)')'Deleted from block ',i , sum(blocks(i)%deleted(:))
          write(*,'(A,x,i2,x,A,x,i8)')&
&        '# of atoms in Block',i,':',blocks(i)%n_atoms-sum(blocks(i)%deleted(:))
        maxtype=maxval(blocks(i)%atype(:))
        print*,maxtype,maxtypeo
      end do
      maxtypeo=max(maxtype,maxtypeo)
    end if
    n_nanopartatoms=0
    if( present(particle)) n_nanopartatoms=particle%n_atoms
    n_ceramic_atoms=n_ceramic_atoms-temp
    n_polymer_atoms=0
    do i=1,n_chains
        n_polymer_atoms=n_polymer_atoms+chains(i)%n_mers*base_mers(chains(i)%mer_id)%n_mer_atoms
    end do
    !n_polymer_atoms=sum(chains(:)%n_mers)
    temp=0
    do i=1,n_chains
      radical=base_mers(chains(i)%mer_id)%radical
      do j=1,chains(i)%n_mers
        do q=1,chains(i)%mers(j)%n_mer_atoms
          temp=temp+chains(i)%mers(j)%deleted(q)
        end do
      end do
      if ( chains(i)%chaintype .eq. 'amorph'.and. .not. radical) temp=temp-2 ! add tail/head H atoms
    end do
    n_polymer_atoms=n_polymer_atoms-temp ! subtract deleted!
    n_total_atoms=n_polymer_atoms+n_ceramic_atoms+n_nanopartatoms
    maxtype=1
    print*, '---'
    allocate(bondtable(n_polymer_atoms*5,2))
    j=1
    do i=1,size(base_mers)
        j=maxval(base_mers(i)%atype(:))
        if ( j > maxtype ) maxtype=j
    end do
    maxtypeo=max(maxtype,maxtypeo)
    allocate(printedatoms(n_total_atoms))
    write(6,*) maxtypeo
    if( n_polymer_atoms  > 0 ) &
    write(*,'(A,x,i6)')'Polymer Atoms:',n_polymer_atoms
    if( present(particle))&
&    write(*,'(A,x,i8)')'Nanoparticle Atoms :',n_nanopartatoms
    write(*,'(A,x,i8)')'Total Number of atoms : ',n_total_atoms
    open(144,file=filelammps)
    open(145,file=filecfg)
    open(146,file=filexyz)
    open(147,file=filebgf)
    write(144,1001)n_polymer_atoms,n_ceramic_atoms
    write(144,*)
    write(144,1002) n_total_atoms
    write(144,1003)
    write(144,1004) maxtypeo
    write(144,*)
    write(144,1005) 0.d0 , BOX(1)
    write(144,1006) 0.d0 , BOX(2)
    write(144,1007) 0.d0 , BOX(3)
    write(144,*)
    write(144,1008)
    write(144,*)
    write(145,'(A,x,i6)') 'Number of particles = ',n_total_atoms
    write(145,'(A)') 'A =1.0'
    write(145,'(A,x,f13.6)')'H0(1,1)=', BOX(1)
    write(145,'(A,x,f13.6)')'H0(1,2)=', 0.0 
    write(145,'(A,x,f13.6)')'H0(1,3)=', 0.0 
    write(145,'(A,x,f13.6)')'H0(2,1)=', 0.0 
    write(145,'(A,x,f13.6)')'H0(2,2)=', BOX(2)
    write(145,'(A,x,f13.6)')'H0(2,3)=', 0.0 
    write(145,'(A,x,f13.6)')'H0(3,1)=', 0.0 
    write(145,'(A,x,f13.6)')'H0(3,2)=', 0.0 
    write(145,'(A,x,f13.6)')'H0(3,3)=', BOX(3)
    write(145,'(A)') '.NO_VELOCITY.'
    write(145,'(A)')'entry_count = 4 '
    write(145,'(A)')'auxiliary[0] = q'
    write(146,*) n_total_atoms
    write(146,*)
    write(147,'(A)')'XTLGRF 200'
    write(147,'(A,x,A)')'DESCRP',filename
    write(147,'(A)')'REMARK BGF FILE CREATED WITH POLYBUILD'
    write(147,'(A,2x,6(f10.5,x))') 'CRYSTX',BOX(1:3),90.0,90.0,90.0
    write(147,'(A)')'FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)'
205 format(a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)        
    cnt=0
    cbond=0
    do i=1,n_boxes
      do j=1,blocks(i)%n_atoms
        if( blocks(i)%deleted(j) .eq. 0) then
          cnt=cnt+1
          write(144,1009) cnt,blocks(i)%atype(j),0.d0, blocks(i)%coord(j,:)
          write(145,'(f12.4)') blocks(i)%mass(j)
          write(145,'(A3)') blocks(i)%symbol(j)
          rvec(1)=blocks(i)%coord(j,1)/BOX(1)
          rvec(2)=blocks(i)%coord(j,2)/BOX(2)
          rvec(3)=blocks(i)%coord(j,3)/BOX(3)
          write(145,'(4(f13.6,x))') rvec,1.d0
          write(146,*) blocks(i)%symbol(j),blocks(i)%coord(j,:)
          write(147,205)'HETATM',cnt,blocks(i)%symbol(j),'','','',&
&              blocks(i)%coord(j,:),blocks(i)%symbol(j),0,0,0.0
          printedatoms(cnt)%symbol=blocks(i)%symbol(j)
          printedatoms(cnt)%charge=blocks(i)%charge(j)
          printedatoms(cnt)%atomtype=blocks(i)%atype(j)
          printedatoms(cnt)%coord=blocks(i)%coord(j,:)
          printedatoms(cnt)%molid=i
          printedatoms(cnt)%index=cnt

        end if
      end do
    end do
    do i=1,n_chains+nchains_functional
      do j=1,chains(i)%n_mers
        do k=1,chains(i)%mers(j)%n_mer_atoms
          if(chains(i)%mers(j)%deleted(k) .eq. 0) then
            cnt=cnt+1
            chains(i)%mers(j)%atom_index(k)=cnt
            write(144,1009) cnt,chains(i)%mers(j)%atype(k),&
&                         0.d0, chains(i)%mers(j)%coord(k,:)
            write(145,'(f12.4)') chains(i)%mers(j)%mass(k)
            write(145,'(A3)') chains(i)%mers(j)%symbol(k)
            rvec(1)=chains(i)%mers(j)%coord(k,1)/BOX(1)
            rvec(2)=chains(i)%mers(j)%coord(k,2)/BOX(2)
            rvec(3)=chains(i)%mers(j)%coord(k,3)/BOX(3)
            write(145,'(4(f13.6,x))') rvec,1.d0
            write(146,*) chains(i)%mers(j)%symbol(k),chains(i)%mers(j)%coord(k,:)
            write(147,205)'HETATM',cnt,chains(i)%mers(j)%symbol(k),&
&                         '','','',chains(i)%mers(j)%coord(k,:),&
&                                chains(i)%mers(j)%symbol(k),0,0,0.0
            printedatoms(cnt)%symbol=chains(i)%mers(j)%symbol(k)
            printedatoms(cnt)%atomtype=chains(i)%mers(j)%atype(k)
            printedatoms(cnt)%charge=chains(i)%mers(j)%charge(k)
            printedatoms(cnt)%coord=chains(i)%mers(j)%coord(k,:)
            printedatoms(cnt)%molid=n_boxes+i
            printedatoms(cnt)%index=cnt
            printedatoms(cnt)%mer_index = j
            printedatoms(cnt)%atom_mer_index = k
          end if
        end do
      end do
                 

               
      radical=base_mers(chains(i)%mer_id)%radical
      if( chains(i)%chaintype .eq. 'amorph' .and. .not. radical) then
        cnt=cnt+1
        write(144,1009) cnt,base_mers(chains(i)%mer_id)%htype, 0.d0, chains(i)%head(:)
        rvec=chains(i)%head/box
        write(145,'(f12.4)') 1.0
        write(145,'(A3)') 'H '
        write(145,'(4(f13.6,x))') rvec,1.d0
        write(146,*) 'H  ',chains(i)%head(:)
        write(147,205)'HETATM',cnt,'H  ','','','',chains(i)%head(:),'H  ',0,0,0.0
            printedatoms(cnt)%symbol='H ' 
            printedatoms(cnt)%atomtype=base_mers(chains(i)%mer_id)%htype
            printedatoms(cnt)%charge=base_mers(chains(i)%mer_id)%headcharge
            printedatoms(cnt)%coord=chains(i)%head(:)
            printedatoms(cnt)%molid=n_boxes+i
            printedatoms(cnt)%index=cnt
        chains(i)%headh_atom_index = cnt
        cnt=cnt+1
        write(144,1009) cnt,base_mers(chains(i)%mer_id)%htype, 0.d0,chains(i)%tail(:)
        rvec=chains(i)%tail/box
        write(145,'(f12.4)') 1.0
        write(145,'(A3)') 'H '
        write(145,'(4(f13.6,x))') rvec,1.d0
        write(146,*) 'H  ',chains(i)%tail(:)
        write(147,205)'HETATM',cnt,'H  ','','','',chains(i)%tail(:),'H  ',0,0,0.0

            printedatoms(cnt)%symbol='H ' 
            printedatoms(cnt)%atomtype=base_mers(chains(i)%mer_id)%htype
            printedatoms(cnt)%charge=base_mers(chains(i)%mer_id)%tailcharge
            printedatoms(cnt)%coord=chains(i)%tail(:)
            printedatoms(cnt)%molid=n_boxes+i
            printedatoms(cnt)%index=cnt
        chains(i)%tailh_atom_index = cnt
        
      end if

      !----------------------------------------------------
      ! Connection table for chain i
      ! ---------------------------------------------------
      do j=1, chains(i)%n_mers
         do k=1,chains(i)%mers(j)%n_mer_atoms
            iatom=chains(i)%mers(j)%atom_index(k)
            do b=1, chains(i)%mers(j)%bondlist(k,1)
               q= chains(i)%mers(j)%bondlist(k,1+b)
               jatom=chains(i)%mers(j)%atom_index(q)
               if( iatom > jatom) then
                   cbond=cbond+1
                   bondtable(cbond,1)=iatom
                   bondtable(cbond,2)=jatom 
                   !print*,cbond,bondtable(cbond,:)
               end if
            end do
         end do
         if( j < chains(i)%n_mers) then
             tc= chains(i)%mers(j)%tailc
             hc= chains(i)%mers(j+1)%headc
             iatom =chains(i)%mers(j)%atom_index(tc)
             jatom =chains(i)%mers(j+1)%atom_index(hc)
             cbond=cbond+1
             bondtable(cbond,1)=iatom
             bondtable(cbond,2)=jatom
             !print*,cbond,bondtable(cbond,:)
         end if
      end do
      ! connect head and tail hydrogens
      if( chains(i)%chaintype .eq. 'amorph' .and. .not. radical) then
          iatom = chains(i)%headh_atom_index
          j = chains(i)%mers(1)%headc
          jatom = chains(i)%mers(1)%atom_index(j)
          cbond=cbond+1
          bondtable(cbond,1) = iatom
          bondtable(cbond,2) =jatom
          !print*,cbond,bondtable(cbond,:)

          iatom = chains(i)%tailh_atom_index
          j = chains(i)%mers(chains(i)%n_mers)%tailc
          jatom = chains(i)%mers(chains(i)%n_mers)%atom_index(j)
          cbond=cbond+1
          bondtable(cbond,1) = iatom
          bondtable(cbond,2) =jatom
         !print*,cbond,bondtable(cbond,:)
     end if 




    end do
    n_total_bonds= cbond
    if( present(particle))then
      cnt=cnt+1
      do i=1,n_nanopartatoms
        write(144,1009) cnt,particle%atype(i),0.d0 , particle%coord(i,:)
        write(145,'(f12.4)') particle%mass(i)
        write(145,'(A3)') particle%symbol(i)
        rvec(1)=particle%coord(i,1)/BOX(1)
        rvec(2)=particle%coord(i,2)/BOX(2)
        rvec(3)=particle%coord(i,3)/BOX(3)
        write(145,'(4(f13.6,x))') rvec,1.d0
        write(146,*) particle%symbol(i), particle%coord(i,:)
        write(147,205)'HETATM',cnt,particle%symbol(i),'','','',&
&                   particle%coord(i,:),particle%symbol(i),0,0,0.0
      end do
    end if      
    iostat=0
        !-----------------------------------------------------------
    close(145)
    close(146)
    close(147)
1001    format('Generated with NewCompositeBuilder ',i10,' Polymer Atoms ',i10,&
            &           ' Ceramic Atoms ')
1002    format(i10,' atoms ')
1003    format(5x,'0 impropers')
1004    format(i3,5x,' atom types')
1005    format(2(f16.4,2x) , 'xlo xhi')
1006    format(2(f16.4,2x) , 'ylo yhi')
1007    format(2(f16.4,2x) , 'zlo zhi')
1008    format('Atoms')
1009    format(i10,2x,i2,3x,f8.2,3f16.6)
  end function print_all
    !----------------------------
  function print_chains(chains,filename) result(iostat)
    type(chain) :: chains(:)
    character(len=*):: filename
    integer :: iostat
    integer :: i, j,k,cnt,n_polymer_atoms
    real*8 , dimension(3) :: vec
    n_polymer_atoms=0  !sum(chains(:)%n_mers)
    do i=1,size(chains)
        n_polymer_atoms=n_polymer_atoms+base_mers(chains(i)%mer_id)%n_mer_atoms*chains(i)%n_mers
    end do
      n_polymer_atoms= n_polymer_atoms+n_amorph_chains*2
      
    cnt=0
    open(144,file=filename)
    write(144,1001)n_polymer_atoms
    write(144,*)
    write(144,1002) n_polymer_atoms
    write(144,*)
    write(144,1004)
    write(144,*)
    write(144,1005) 0.d0 , BOX(1)
    write(144,1006) 0.d0 , BOX(2)
    write(144,1007) 0.d0 , BOX(3)
    write(144,*)
    write(144,1008)
    write(144,*)
    do i=1,size(chains(:))
      !print*,' chain no : ', i
      if( chains(i)%chaintype .eq. 'amorph') then
        vec=chains(i)%coord(1,:)-chains(i)%coord(2,:)
        vec=vec/sqrt(dot_product(vec,vec))
        vec=1.15d0*vec
        vec=chains(i)%coord(1,:)+vec
        vec=chains(i)%head(:)
        ! vec is the position of chain head cap hydrogen!
        cnt=cnt+1
        ! Hydrogen type for reax is 2 !
        write(144,1009)cnt,2 , 0.d0   ,vec
      end if
      do j=1,chains(i)%n_mers
        do k=1,chains(i)%mers(j)%n_mer_atoms
          cnt=cnt+1
          write(144,1009) cnt,chains(i)%mers(j)%atype(k),&
&                       0.d0, chains(i)%mers(j)%coord(k,:)
        end do
      end do
      if( chains(i)%chaintype .eq. 'amorph') then
        vec=chains(i)%coord(chains(i)%n_mers,:)-chains(i)%coord(chains(i)%n_mers-1,:)
        vec=vec/sqrt(dot_product(vec,vec))
        vec=1.15d0*vec
        vec=chains(i)%coord(chains(i)%n_mers,:)+vec
      ! vec is the position of chain head cap hydrogen!
        vec=chains(i)%tail(:)
        cnt=cnt+1
      ! Hydrogen type for reax is 2 !
        write(144,1009)cnt,2 , 0.d0   ,chains(i)%tail(:)
      end if
    end do
    iostat=0
1001    format('Generated with NewCompositeBuilder ',i10,' Polymer Atoms ',i10,&
            &           ' Ceramic Atoms ')
1002    format(i10,' atoms ')
1004    format(5x,'4 atom types')
1005    format(2(f16.4,2x) , 'xlo xhi')
1006    format(2(f16.4,2x) , 'ylo yhi')
1007    format(2(f16.4,2x) , 'zlo zhi')
1008    format('Atoms')
1009    format(i10,2x,i2,3x,f8.2,3f16.6)
  end function print_chains
!===================================================================================
  function symbol_to_hell(symbol) result(hellid)
        character(len=*) :: symbol
        integer :: hellid

        hellid=0
        if( symbol .eq. 'O ') hellid = 1
        if( symbol .eq. 'C ') hellid = 2
        if( symbol .eq. 'H ') hellid = 3
        if( symbol .eq. 'Al') hellid = 4
    end function symbol_to_hell



end module print_files
