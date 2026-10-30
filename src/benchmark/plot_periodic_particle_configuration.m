function plot_periodic_particle_configuration(ax,T,Lx,Ly)

    % Draw actual particle centers/radii in the primary periodic cell.
    % Primary particles are plotted in blue; periodic image portions that
    % cross an edge are shown in light gray outside/through the opposite edge.

    x = T.x(:);
    y = T.y(:);
    r = T.radius(:);

    primaryFace = [0.82 0.88 0.96];
    primaryEdge = [0.10 0.35 0.70];
    imageFace   = [0.92 0.92 0.92];
    imageEdge   = [0.55 0.55 0.55];

    % Periodic image copies first.
    shifts = [-1 0 1];
    for i = 1:numel(r)
        for sx = shifts
            for sy = shifts
                if sx == 0 && sy == 0
                    continue;
                end

                xc = x(i) + sx*Lx;
                yc = y(i) + sy*Ly;
                rr = r(i);

                if (xc+rr >= 0) && (xc-rr <= Lx) && ...
                        (yc+rr >= 0) && (yc-rr <= Ly)
                    rectangle(ax,'Position',[xc-rr,yc-rr,2*rr,2*rr], ...
                        'Curvature',[1 1], ...
                        'FaceColor',imageFace, ...
                        'EdgeColor',imageEdge, ...
                        'LineWidth',0.35);
                end
            end
        end
    end

    % Primary particles.
    for i = 1:numel(r)
        rectangle(ax,'Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)], ...
            'Curvature',[1 1], ...
            'FaceColor',primaryFace, ...
            'EdgeColor',primaryEdge, ...
            'LineWidth',0.55);
    end

    % Primary periodic cell boundary.
    rectangle(ax,'Position',[0 0 Lx Ly], ...
        'EdgeColor','k','LineWidth',1.4);

    axis(ax,'equal');
    xlim(ax,[0 Lx]);
    ylim(ax,[0 Ly]);
    set(ax,'XTick',[],'YTick',[]);
    xlabel(ax,'');
    ylabel(ax,'');
end
